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

#include "EMSolver.h"
#include "Helpers.h"

EMSolver::EMSolver()
{
	ga = NULL;
	iterative = NULL;
	status = Clean;
}

EMSolver::~EMSolver()
{
	delete ga;
	delete iterative;
}

bool EMSolver::Formulate(const DataStore &store)
{
	if (status != Clean)
		return false;
	ContigCount = store.ContigCount;
	U.resize(ContigCount, true);
	T.resize(ContigCount);
	X.resize(ContigCount);
	this->store = store;
	status = Formulated;
	return true;
}

bool EMSolver::Solve()
{
	if (status < Formulated)
		return false;
	timerId = Helpers::ElapsedTimers.AddTimer();
	bestObjective = -Helpers::Inf;
	Iteration = 0;
	prepareSlacks();
	bool success = true;
	double lastElapsed = 0;
	while (!shouldTerminate())
	{
		if (Iteration > 0)
			lastT = ga->T;
		if (!expectation())
		{
			if (Options.VerboseOutput > 0)
				fprintf(stderr, "   [+] Iteration %2i: failed expectation (%9.2lf ms).\n", Iteration, getTime(lastElapsed));
			success = false;
			break;
		}
		if (Options.VerboseOutput > 0)
			fprintf(stderr, "   [+] Iteration %2i: expectation  of %10.3lf in %9.2lf ms.\n", Iteration, ga->GetObjective(), getTime(lastElapsed));
		if (!maximization())
		{
			if (Options.VerboseOutput > 0)
				fprintf(stderr, "   [+] Iteration %2i: failed maximization (%9.2lf ms).\n", Iteration, getTime(lastElapsed));
			success = false;
			break;
		}
		if (Options.VerboseOutput > 0)
			fprintf(stderr, "   [+] Iteration %2i: maximization of %10.3lf in %9.2lf ms.\n", Iteration, iterative->GetObjective(), getTime(lastElapsed));
		fprintf(stderr, "\n");
		updateSlack();
		updateBest();
		Iteration++;
	}
	if (success)
	{
		T = bestT;
		X = bestX;
		objectiveValue = bestObjective;
	}
	status = (success ? Success : Fail);
	Helpers::ElapsedTimers.RemoveTimer(timerId);
	return success;
}

SolverStatus EMSolver::GetStatus() const
{
	return status;
}

double EMSolver::GetObjective() const
{
	if (status == Success)
		return objectiveValue;
	return -Helpers::Inf;
}

IloAlgorithm::Status EMSolver::GetCplexStatus() const
{
	if (iterative != NULL)
		return iterative->GetCplexStatus();
	return (IloAlgorithm::Status)0;
}

double EMSolver::GetDistanceSlack(int i) const
{
	if (iterative != NULL)
		return iterative->GetDistanceSlack(i);
	return -1;
}

double EMSolver::GetOrderSlack(int i) const
{
	if (iterative != NULL)
		return iterative->GetOrderSlack(i);
	return -1;
}

int EMSolver::GetSlackCount() const
{
	if (iterative != NULL)
		return iterative->GetSlackCount();
	return -1;
}

const DataStore &EMSolver::GetStore() const
{
	return store;
}

void EMSolver::prepareSlacks()
{
	int count = 0;
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
		count++;
	distanceSlack.assign(count, 0);
	orderSlack.assign(count, 0);
}

bool EMSolver::converged()
{
	return Iteration > 0 && lastT == ga->T;
}

bool EMSolver::shouldTerminate()
{
	if (converged())
		return true;
	if (Options.TimeLimit > 0 && Helpers::ElapsedTimers.Elapsed(timerId) > (double)Options.TimeLimit * 1000)
		return true;
	return false;
}

bool EMSolver::expectation()
{
	delete ga;
	ga = new GASolver();
	ga->Options = Options;
	if (!ga->Formulate(store, distanceSlack, orderSlack))
		return false;
	if (iterative != NULL && iterative->GetStatus() == Success)
		ga->AddIndividual(iterative->T);
	if (!ga->Solve())
		return false;
	return true;
}

bool EMSolver::maximization()
{
	delete iterative;
	iterative = new IterativeSolver(ga->U, ga->T, ContigCount);
	iterative->Options = Options;
	if (!iterative->Formulate(store))
		return false;
	if (!iterative->Solve())
		return false;
	return true;
}

void EMSolver::updateSlack()
{
	int num = 0, id = 0;
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++, id++)
	{
		int a = it->first.first, b = it->first.second;
		if ((ga->T[a] ^ ga->T[b]) != it->second.EqualOrientation)
		{
			distanceSlack[id] = iterative->GetDistanceSlack(num);
			orderSlack[id] = iterative->GetOrderSlack(num);
			num++;
		}
	}
}

void EMSolver::updateBest()
{
	double newObjective = iterative->GetObjective();
	if (newObjective > bestObjective)
	{
		bestObjective = newObjective;
		bestT = iterative->T;
		bestX = iterative->X;
	}
}

double EMSolver::getTime(double &lastIteration) const
{
	double last = lastIteration;
	lastIteration = Helpers::ElapsedTimers.Elapsed(timerId);
	return lastIteration - last;
}
