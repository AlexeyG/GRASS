#ifndef _FIXEDMIQPSOLVER_H
#define _FIXEDMIQPSOLVER_H

#include "Solver.h"
#include "DataStore.h"
#include "SolverConfiguration.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

class FixedMIQPSolver : public Solver
{
public:
	FixedMIQPSolver(const vector<bool> &u, const vector<bool> &t, int len);
	virtual ~FixedMIQPSolver();

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

private:
	bool addContigs(const DataStore &store);
	bool addContig(const Contig &contig);
	bool addLinks(const DataStore &store);
	bool addLink(int a, int b, const ContigLink &link);
	bool addDistanceConstraint(int a, int b, bool e, bool r, double sigma, double mu, IloNumVar &xi_l);
	bool addOrderConstraint(int a, int b, bool e, bool r, IloNumVar &delta_l);
	void appendOrientationObjective(int a, int b, bool e, double w);
	bool appendDistanceObjective(int a, int b, bool e, double w, const IloNumVar &xi_l);
	bool appendOrderObjective(int a, int b, bool e, double w, const IloNumVar &delta_l);
	void appendSizeObjective();
	bool createModel();
	void saveSolution();

public:
	SolverConfiguration CPLEXOptions;

public:
	static const double SlackMax = 5e5;
	static const double CoordMax = 1e10;
	static const double DesiredDistanceSlackMax = 6;
	static const double DesiredOrderSlackMax = 6;

private:
	IloEnv environment;
	IloModel model;
	IloNumVarArray x;
	IloNumVarArray xi, delta;
	IloRangeArray constraints;
	IloExpr h,p;
	double g,s;
	IloCplex cplex;
	double bestObjective;
	vector<int> len;
	vector<bool> optimized;
};

#endif
