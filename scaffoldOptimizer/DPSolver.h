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


#ifndef _DPSOLVER_H
#define _DPSOLVER_H
#include "Solver.h"
#include "DPGraph.h"
#include "SolverConfiguration.h"
#include "ScaffoldExtractor.h"

#include <vector>

using namespace std;

class DPSolver : public Solver
{
public:
	DPSolver();
	virtual ~DPSolver();

public:
	virtual bool Formulate(const DataStore &store);
	virtual bool Solve();
	virtual SolverStatus GetStatus() const;
	virtual double GetObjective() const;

public:
	vector<Scaffold> Scaffolds;

public:
	const static int ScaffoldSeprator = 10;

private:
	bool processComponents();

public:
	int MaxIteration;

private:
	DataStore store;
	DPGraph graph;
	vector< vector<int> > connectedComponents;
	int nComponents;
	double objectiveValue;
};
#endif
