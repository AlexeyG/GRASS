#ifndef _EXPECTATIONMAXIMIZATION_H
#define _EXPECTATIONMAXIMIZATION_H
#include "Solver.h"
#include "GASolver.h"

class ExpectationMaximization : Solver
{
public:
	ExpectationMaximization();
	virtual ~ExpectationMaximization();

public:
	bool Formulate(const DataStore &store, const vector<double> &coord);
	virtual bool Formulate(const DataStore &store);
	virtual bool Solve();
	virtual SolverStatus GetStatus() const;
	virtual double GetObjective() const;

public:
	double GetDistanceSlack(int i) const;
	double GetOrderSlack(int i) const;
	int GetSlackCount() const;

protected:
	SolverStatus status;
};
#endif
