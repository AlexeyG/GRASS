#ifndef _BRANCHANDBOUND_H
#define _BRANCHANDBOUND_H

#include "Solver.h"
#include "DataStore.h"
#include "SolverConfiguration.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

class BranchAndBound : public Solver
{
public:
	BranchAndBound(const vector<bool> &u, const vector<bool> &t, int len);
	virtual ~BranchAndBound();

public:
	bool Formulate(const DataStore &store, const vector<double> &coord);
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
	bool formulate(const DataStore &store);
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
	bool addCoordinateConstraints(const vector<double> &coord);
	bool createModel();
	bool assignPriorities(const DataStore &store);
	void saveSolution();

public:
	SolverConfiguration CPLEXOptions;

public:
	static const double SlackMax = 5e7;
	static const double CoordMax = 1e10;
	//static const double DesiredDistanceSlackMax = 6;
	//static const double DesiredOrderSlackMax = 1000;
	static const double DesiredDistanceSlackMax = 6;
	static const double DesiredOrderSlackMax = 1;

public:
	vector<bool> Incumbent;
	vector<double> Slack;

private:
public:
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
