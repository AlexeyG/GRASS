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

#include "DPSolver.h"
#include "EMSolver.h"
#include "Helpers.h"
#include "MinMax.h"

DPSolver::DPSolver()
	: nComponents(0)
{
	ContigCount = 0;
	status = Clean;
	objectiveValue = 0;
	MaxIteration = 0;
}

DPSolver::~DPSolver()
{
}

bool DPSolver::Formulate(const DataStore &store)
{
	if (status != Clean)
		return false;
	ContigCount = store.ContigCount;
	U.resize(ContigCount);
	T.resize(ContigCount);
	X.resize(ContigCount);
	this->store = store;
	graph = DPGraph(store);
	nComponents = graph.FindConnectedComponents(connectedComponents);
	status = Formulated;
	return true;
}

bool DPSolver::Solve()
{
	bool result;
	if (status < Formulated)
		return false;
	fprintf(stderr, "[i] Have %i connected components.\n", nComponents);
	if (nComponents == 0)
		result = true;
	else
		result = processComponents();
	status = (result ? Success : Fail);
	return result;
}

SolverStatus DPSolver::GetStatus() const
{
	return status;
}

double DPSolver::GetObjective() const
{
	if (status == Success)
		return objectiveValue;
	return -Helpers::Inf;
}

bool DPSolver::processComponents()
{
	bool result = true;
	vector<double> minX(nComponents);
	vector<double> maxX(nComponents);
	vector< vector<int> > backTransform(nComponents);
	vector<Scaffold> scaffolds;
	MaxIteration = 0;
	for (int i = 0; i < nComponents; i++)
	{
		int nContigsComponent = connectedComponents[i].size();
		fprintf(stderr, "    [i] Processing component %i of size %i.\n", i + 1, nContigsComponent);
		DataStore compStore;
		EMSolver *solver = new EMSolver();
		solver->Options = Options;
		store.Extract(connectedComponents[i], compStore, backTransform[i]);
		if (!solver->Formulate(compStore) || !solver->Solve())
		{
			fprintf(stderr, "        [-] Unable to solve or formulate.\n");
			result = false;
			delete solver;
			break;
		}
		else
			fprintf(stderr, "        [+] Formulated and solved subproblem.\n");
		scaffolds = ScaffoldExtractor::Extract(*solver);
		for (vector<Scaffold>::iterator it = scaffolds.begin(); it != scaffolds.end(); it++)
		{
			it->ApplyTransform(backTransform[i]);
			it->NormalizeCoordindates();
		}
		Scaffolds.insert(Scaffolds.end(), scaffolds.begin(), scaffolds.end());
		objectiveValue += solver->GetObjective();
		MaxIteration = max(MaxIteration, solver->Iteration);
		minX[i] =   Helpers::Inf;
		maxX[i] = - Helpers::Inf;
		for (int j = 0; j < nContigsComponent; j++)
		{
			int id = backTransform[i][j];
			U[id] = solver->U[j];
			T[id] = solver->T[j];
			X[id] = solver->X[j];
			if (solver->U[j])
			{
				int contigLen = compStore[j].Sequence.Nucleotides.length();
				minX[i] = min(minX[i], (solver->T[j] == 1 ? solver->X[j] - contigLen + 1 : solver->X[j]));
				maxX[i] = max(maxX[i], (solver->T[j] == 0 ? solver->X[j] + contigLen - 1 : solver->X[j]));
			}
		}
		delete solver;
	}
        
        if (result)
        {
            fprintf(stderr, "    [i] Adjusting subscaffold positions\n");
            for (int i = 0; i < nComponents; i++)
            {
                    int shift = minX[i];
                    int offset = (i > 0 ? maxX[i - 1] + ScaffoldSeprator : 0); // that's a strange statement - was I planning to put all scaffolds on a single line?
                    //int offset = 0;
                    int nContigsComponent = connectedComponents[i].size();
                    for (int j = 0; j < nContigsComponent; j++)
                            if (U[backTransform[i][j]])
                                    X[backTransform[i][j]] -= shift - offset;
            }
        }
	return result;
}
