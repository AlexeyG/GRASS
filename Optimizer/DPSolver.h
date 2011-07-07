#ifndef _DPSOLVER_H
#define _DPSOLVER_H
#include "Solver.h"
#include "DPGraph.h"
#include "SolverConfiguration.h"

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
	const static int ScaffoldSeprator = 10;

private:
	bool processComponents();
	bool processSingleProblem();

private:
	DataStore store;
	DPGraph graph;
	vector< vector<int> > connectedComponents;
	int nComponents;
	double objectiveValue;
};
#endif
