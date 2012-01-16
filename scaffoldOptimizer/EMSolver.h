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


#ifndef _EMSOLVER_H
#define _EMSOLVER_H
#include "Solver.h"
#include "GASolver.h"
#include "IterativeSolver.h"

class EMSolver : public Solver
{
public:
	EMSolver();
	~EMSolver();

public:
	virtual bool Formulate(const DataStore &store);
	virtual bool Solve();
	virtual SolverStatus GetStatus() const;
	virtual double GetObjective() const;

public:
	IloAlgorithm::Status GetCplexStatus() const;
	double GetDistanceSlack(int i) const;
	double GetOrderSlack(int i) const;
	int GetSlackCount() const;
	const DataStore &GetStore() const;

private:
	void prepareSlacks();
	bool converged();
	bool shouldTerminate();
	bool expectation();
	bool maximization();
	void updateSlack();
	void updateBest();
	double getTime(double &lastIteration) const;

protected:
	SolverStatus status;

public:
	int Iteration;

private:
	int timerId;
	GASolver *ga;
	IterativeSolver *iterative;
	DataStore store;
	vector<double> distanceSlack, orderSlack;
	vector<bool> bestT;
	vector<double> bestX;
	vector<bool> lastT;
	double bestObjective;
	double objectiveValue;
};

#endif
