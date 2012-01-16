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

#ifndef _SOLVER_H
#define _SOLVER_H
#include "DataStore.h"
#include "SolverConfiguration.h"
#include <vector>

enum SolverStatus {Clean, Formulated, Success, Fail};

class Solver
{
public:
	Solver() {};
	virtual ~Solver() {};
	virtual bool Formulate(const DataStore &store) = 0;
	virtual bool Solve() = 0;
	virtual SolverStatus GetStatus() const = 0;
	virtual double GetObjective() const = 0;

public:
	vector<double> X;
	vector<bool> U,T;
	int ContigCount;
	SolverConfiguration Options;

protected:
	SolverStatus status;
};
#endif
