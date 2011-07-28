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

private:
	int timerId;
	int iteration;
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
