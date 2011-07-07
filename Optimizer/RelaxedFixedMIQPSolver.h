#ifndef _RELAXEDFIXEDMIQPSOLVER_H
#define _RELAXEDFIXEDMIQPSOLVER_H

#include "Solver.h"
#include "DataStore.h"
#include "SolverConfiguration.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

class RelaxedFixedMIQPSolver : public Solver
{
public:
	RelaxedFixedMIQPSolver(const vector<bool> &u, const vector<bool> &t, int len);
	virtual ~RelaxedFixedMIQPSolver();

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
	bool addLink(int a, int b, const ContigLink &link, int &num);
	bool addDistanceConstraint(int a, int b, bool e, bool r, double sigma, double mu);
	bool addOrderConstraint(int a, int b, bool e, bool r);
	void appendOrientationObjective(int a, int b, bool e, double w);
	bool appendDistanceObjective(int a, int b, bool e, double w, int num);
	bool appendOrderObjective(int a, int b, bool e, double w, int num);
	void appendSizeObjective();
	bool createModel();
	void saveSolution();

public:
	SolverConfiguration CPLEXOptions;

public:
	static const double SlackMax = 5e5;
	static const double CoordMax = 1e10;
	static const double WeightMultiplier = 1e7;
	static const double DesiredDistanceSlackMax = 6;
	static const double DesiredOrderSlackMax = 6;

public:
	vector<bool> Incumbent;
	vector<double> Slack;

private:
	IloEnv environment;
	IloModel model;
	IloNumVarArray x;
	IloNumVarArray xi, delta;
	IloNumVarArray alpha, beta;
	IloNumVarArray xi_alpha, delta_beta;
	IloRangeArray constraints;
	IloExpr h,p;
	double g,s;
	IloCplex cplex;
	vector<int> len;
	vector<bool> optimized;
};

#endif
