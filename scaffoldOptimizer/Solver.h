#ifndef _SOLVER_H
#define _SOLVER_H
#include "DataStore.h"
#include "SolverConfiguration.h"
#include <vector>

enum SolverStatus {Clean, Formulated, Success, Fail};

class Solver
{
public:
	Solver();
	virtual ~Solver();
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
