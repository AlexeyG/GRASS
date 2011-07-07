#ifndef _ITERATIVESOLVER_H
#define _ITERATIVESOLVER_H
#include "Solver.h"
#include "DataStore.h"
#include "FixedMIQPSolver.h"
#include "ExtendedFixedMIQPSolver.h"

class IterativeSolver : public Solver
{
public:
	IterativeSolver(const vector<bool> &u, const vector<bool> &t, int length);
	virtual ~IterativeSolver();
	virtual bool Formulate(const DataStore &store);
	virtual bool Solve();
	virtual SolverStatus GetStatus() const;
	virtual double GetObjective() const;

public:
	IloAlgorithm::Status GetCplexStatus() const;
	double GetDistanceSlack(int i) const;
	double GetOrderSlack(int i) const;
	int GetSlackCount() const;
	double GetHeuristicObjective() const;

public:
	vector<bool> U,T;
	int Disabled;

protected:
	SolverStatus status;

private:
	ExtendedFixedMIQPSolver solver;
	ExtendedFixedMIQPSolver extension;
	DataStore store;
};
#endif
