/*
 * scaffoldOptimizer : solves the MIQP optimization and produces linear scaffold
 * sequences.
 * Copyright (C) 2011  Alexey Gritsenko
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 * 
 * 
 * Email: a.gritsenko@tudelft.nl
 * Mail: Delft University of Technology
 *       Faculty of Electrical Engineering, Mathematics, and Computer Science
 *       Department of Mediamatics
 *       P.O. Box 5031
 *       2600 GA, Delft, The Netherlands
 */

#include "GASolver.h"
#include "Helpers.h"
#include "RandomizedGreedyInitializer.h"
#include "ExtendedFixedMIQPSolver.h"
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <set>
#include "omp.h"

using namespace std;

GASolver::GASolver()
{
	SelectionSize = 40;
	CrossoverRate = 0.5;
	RestartGenerations = 30;
	LocalSearchM = 50;
	status = Clean;
}

GASolver::~GASolver()
{
}

bool GASolver::Formulate(const DataStore &store)
{
	return Formulate(store, vector<double>(), vector<double>());
}

bool GASolver::Formulate(const DataStore &store, const vector<double> &distanceSlack, const vector<double> &orderSlack)
{
	if (status != Clean)
		return false;
	ContigCount = store.ContigCount;
	U.assign(ContigCount, true);
	T.resize(ContigCount);
	X.resize(ContigCount); // compability
	formulateMatrix(store, distanceSlack, orderSlack);
	Options.SuppressOutput = true;
	bestObjective = -Helpers::Inf;
	populationSize = 0;
	status = Formulated;
	return true;
}

bool GASolver::Solve()
{
	if (status < Formulated)
		return false;

	timerId = Helpers::ElapsedTimers.AddTimer();
	iteration = 0;
	restartCount = 0;
	lastSuccess = 0;
	double lastTime = 0;
	
	omp_set_num_threads(Options.Threads);
	localSearch(generatePopulation(populationSize));
	if (Options.VerboseOutput > 1)
		printf("        [i] Generated population: %.2lf ms\n", getTime(lastTime));
	selectInitialSolution();
	while (!shouldTerminate())
	{
		localSearch(crossover());
		select();
		if (Options.VerboseOutput > 1)
			printf("        [i] Iteration %i: %.2lf ms\n", iteration + 1, getTime(lastTime));
		if (iteration - lastSuccess >= RestartGenerations)
		{
			restartCount++;
			lastSuccess = iteration;
			localSearch(restart(1));
			localSearch(generatePopulation(populationSize));
			if (Options.VerboseOutput > 1)
			{
				if (Options.GARestarts > 0)
					printf("            [i] Restarted: attempt %i out of %i\n", restartCount, Options.GARestarts);
				else
					printf("            [i] Restarted: attempt %i\n", restartCount);
			}
		}
		iteration++;
	}
	Helpers::ElapsedTimers.RemoveTimer(timerId);
	if (status != Fail)
		status = Success;
	return status == Success;
}

SolverStatus GASolver::GetStatus() const
{
	return status;
}

double GASolver::GetObjective() const
{
	if (status == Success)
		return bestObjective;
	return -Helpers::Inf;
}

void GASolver::AddIndividual(const vector<bool> &t)
{
	if ((int)population.size() < populationSize + 1)
		population.resize(populationSize + 1);
	population[populationSize++] = GAIndividual(t, matrix);
	updateSolution(population[populationSize - 1]);
}

bool GASolver::shouldTerminate()
{
	if (Options.GATimeLimit > 0 && Helpers::ElapsedTimers.Elapsed(timerId) > (double)Options.GATimeLimit * 1000)
		return true;
	if (Options.GARestarts > 0 && restartCount >= Options.GARestarts)
		return true;
	return false;
}

int GASolver::generatePopulation(int from)
{
	if (populationSize < SelectionSize)
		population.resize(SelectionSize), populationSize = SelectionSize;
	#pragma omp parallel
	{
		srand(time(NULL) ^ omp_get_thread_num());
		#pragma omp for
		for (int i = from; i < SelectionSize; i++)
		{
			RandomizedGreedyInitializer init(ContigCount, matrix);
			population[i] = init.MakeSolution(matrix);
		}
	}
	return from;
}

int GASolver::localSearch(int from)
{
	#pragma omp parallel
	{
		srand(time(NULL) ^ omp_get_thread_num());
		#pragma omp for
		for (int i = from; i < populationSize; i++)
			RandomizedKopt(population[i]);
	}
	return from;
}

int GASolver::crossover()
{
	int count = (int)(CrossoverRate * SelectionSize);
	int newSize = populationSize + count;
	population.resize(newSize);
	#pragma omp parallel
	{
		srand(time(NULL) ^ omp_get_thread_num());
		#pragma omp for
		for (int i = 0; i < count; i++)
		{
			int a = rand() % populationSize, b = rand() % populationSize;
			population[populationSize + i] = InnovativeCrossover(population[a], population[b]);
		}
	}
	populationSize = newSize;
	return newSize - count;
}

void GASolver::select()
{
	int i = 1, j = 1;
	sort(population.begin(), population.end(), greater<GAIndividual>());
	for (; i < populationSize && j < SelectionSize; i++)
	{
		if (population[i] != population[j - 1])
			population[j++] = population[i];
	}
	populationSize = j;
	population.resize(populationSize);
	updateSolution(population[0]);
}

int GASolver::restart(int from)
{
	#pragma omp parallel
	{
		srand(time(NULL) ^ omp_get_thread_num());
		#pragma omp for
		for (int i = from; i < populationSize; i++)
			Mutate(population[i]);
	}
	return from;
}

void GASolver::formulateMatrix(const DataStore &store, const vector<double> &distanceSlack, const vector<double> &orderSlack)
{
	matrix = GAMatrix(ContigCount + 1);
	matrix[ContigCount][ContigCount] = 1;
	int num = 0;
	int distanceCount = distanceSlack.size(), orderCount = orderSlack.size();
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++, num++)
	{
		int i = it->first.first, j = it->first.second;
		double xi = (num < distanceCount ? distanceSlack[num] : 0);
		double delta = (num < orderCount ? orderSlack[num] : 0);
		double w = it->second.Weight;
		double xiP = (xi >= ExtendedFixedMIQPSolver::DesiredDistanceSlackMax ? 1 : xi / ExtendedFixedMIQPSolver::DesiredDistanceSlackMax);
		double deltaP = (delta >= ExtendedFixedMIQPSolver::DesiredOrderSlackMax ? 1 : delta / ExtendedFixedMIQPSolver::DesiredOrderSlackMax);
		w -= (xiP + deltaP) * w / 2;
		if (it->second.EqualOrientation)
		{
			matrix[i][i] -= w;
			matrix[j][j] -= w;
			matrix[i][j] += 2 * w;
			matrix[ContigCount][ContigCount] += w; // constant summand
		}
		else
		{
			matrix[i][i] += w;
			matrix[j][j] += w;
			matrix[i][j] -= 2 * w;
		}
	}
	for (int i = 0; i < ContigCount; i++)
		for (int j = i + 1; j < ContigCount; j++)
			matrix[i][j] = matrix[j][i] = (matrix[i][j] + matrix[j][i]) / 2.0;
	matrix.CalculatePositions();
}

void GASolver::selectInitialSolution()
{
	double best = population[0].GetObjective();
	int bestInd = 0;
	for (int i = 1; i < populationSize; i++)
		if (best < population[i].GetObjective())
		{
			best = population[i].GetObjective();
			bestInd = i;
		}
	updateSolution(population[bestInd]);
}

void GASolver::updateSolution(const GAIndividual &ind)
{
	double objective = ind.GetObjective();
	if (objective > bestObjective + Helpers::Eps)
	{
		if (Options.VerboseOutput > 1)
			cout << "            [+] Best found: " << objective << endl;
		bestObjective = objective;
		for (int i = 0; i < ContigCount; i++)
			T[i] = ind.X[i];
		lastSuccess = iteration;
	}
}

double GASolver::getTime(double &lastIteration) const
{
	double last = lastIteration;
	lastIteration = Helpers::ElapsedTimers.Elapsed(timerId);
	return lastIteration - last;
}

GAIndividual GASolver::InnovativeCrossover(const GAIndividual &p1, const GAIndividual &p2)
{
	GAIndividual offspring(p1);
	int eqCnt = 0, neqCnt = 0;
	vector<int> eq, neq;
	for (int i = 0; i < ContigCount; i++)
		if (p1.X[i] == p2.X[i])
			eq.push_back(i), eqCnt++;
		else
			neq.push_back(i), neqCnt++;
	for (int i = neqCnt; i > 0; i--)
	{
		random_shuffle(neq.begin(), neq.end());
		int p = -1;
		for (int j = 0; j < neqCnt; j++)
			if (offspring.Gain[neq[j]] > Helpers::Eps)
			{
				p = j;
				break;
			}
		if (p >= 0)
		{
			offspring.Flip(neq[p], matrix);
			swap(neq[p], neq[neqCnt - 1]), neqCnt--, neq.resize(neqCnt);
		}
		if (eqCnt > 0)
		{
			int p = -1;
			for (int j = 0; j < eqCnt; j++)
				if (p < 0 || offspring.Gain[eq[p]] < offspring.Gain[eq[j]])
					p = j;
			offspring.Flip(eq[p], matrix);
			swap(eq[p], eq[eqCnt - 1]), eqCnt--, eq.resize(eqCnt);
		}
	}
	return offspring;
}

void GASolver::RandomizedKopt(GAIndividual &ind)
{
	vector<int> perm(ContigCount);
	for (int i = 0; i < ContigCount; i++)
		perm[i] = i;
	vector<bool> used(ContigCount);
	while (true)
	{
		GAIndividual xprev(ind), xbest(ind);
		used.assign(ContigCount, false);
		double gBest = 0, g = 0;
		int unused = ContigCount, lastBest = 0;
		do
		{
			lastBest++;
			random_shuffle(perm.begin(), perm.end());
			for (int i = 0; i < ContigCount; i++)
				if (ind.Gain[perm[i]] > Helpers::Eps)
				{
					g += ind.Gain[perm[i]];
					ind.Flip(perm[i], matrix);
					used[perm[i]] = true;
					unused--;
					if (g > gBest)
					{
						gBest = g;
						xbest = ind;
						lastBest = 0;
					}
				}
			if (unused > 0)
			{
				int p = -1;
				for (int i = 0; i < ContigCount; i++)
					if (!used[i] && (p < 0 || ind.Gain[p] < ind.Gain[i]))
						p = i;
				g += ind.Gain[p];
				ind.Flip(p, matrix);
				used[p] = true;
				unused--;
				if (g > gBest)
				{
					gBest = g;
					xbest = ind;
					lastBest = 0;
				}
			}
		} while (unused > 0 && lastBest >= LocalSearchM);
		if (gBest > Helpers::Eps)
			ind = xbest;
		else
		{
			ind = xprev;
			break;
		}
	}
}

void GASolver::Mutate(GAIndividual &ind)
{
	int vars = ContigCount / 3;
	vector<int> perm(ContigCount);
	for (int i = 0; i < ContigCount; i++)
		perm[i] = i;
	random_shuffle(perm.begin(), perm.end());
	for (int i = 0; i < vars; i++)
		ind.Flip(perm[i], matrix);
}

/*bool GASolver::checkObjective(GAIndividual &ind)
{
	vector<bool> t(ContigCount);
	for (int i = 0; i < ContigCount; i++)
		t[i] = ind.X[i];
	GAIndividual check(t, matrix);
	if (fabs(check.GetObjective() - ind.GetObjective()) > Helpers::Eps)
		return false;
	return true;
}

bool GASolver::checkGains(GAIndividual &ind)
{
	vector<bool> t(ContigCount);
	for (int j = 0; j < ContigCount; j++)
		t[j] = ind.X[j];
	GAIndividual check(t, matrix);
	double objective = check.GetObjective();
	for (int i = 0; i < ContigCount; i++)
	{
		for (int j = 0; j < ContigCount; j++)
			t[j] = ind.X[j];
		t[i] = !t[i];
		GAIndividual inc(t, matrix);
		double gain = inc.GetObjective() - objective;
		if (fabs(gain - ind.Gain[i]) > Helpers::Eps)
		{
			cout << "Error at " << i << endl;
			cout << ind.Gain[i] << " vs. " << gain << endl;
			return false;
		}
	}
	return true;
}*/
