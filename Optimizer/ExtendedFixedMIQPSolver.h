#ifndef _EXTENDEDFIXEDMIQPSOLVER_H
#define _EXTENDEDFIXEDMIQPSOLVER_H

#include "Solver.h"
#include "DataStore.h"
#include "SolverConfiguration.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

class ExtendedFixedMIQPSolver : public Solver
{
public:
	ExtendedFixedMIQPSolver(const vector<bool> &u, const vector<bool> &t, int len);
	virtual ~ExtendedFixedMIQPSolver();

public:
	bool Formulate(const DataStore &store, const vector<bool> &enabledDistance, const vector<bool> &enabledOrder);
	bool Formulate(const DataStore &store, const vector<bool> &enabledDistance, const vector<bool> &enabledOrder, const vector<double> &coord);
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
	bool formulate(const DataStore &store, const vector<bool> &enabledDistance, const vector<bool> &enabledOrder);
	bool addContigs(const DataStore &store);
	bool addContig(const Contig &contig);
	bool addLinks(const DataStore &store, const vector<bool> &enabledDistance, const vector<bool> &enabledOrder);
	bool addLink(int a, int b, const ContigLink &link, int &num, bool enabledDistance, bool enabledOrder);
	bool addDistanceConstraint(int a, int b, bool e, bool r, double sigma, double mu, IloNumVar &xi_l);
	bool addOrderConstraint(int a, int b, bool e, bool r, IloNumVar &delta_l);
	void appendOrientationObjective(int a, int b, bool e, double w);
	bool appendDistanceObjective(int a, int b, bool e, double w, const IloNumVar &xi_l, bool enabled);
	bool appendOrderObjective(int a, int b, bool e, double w, const IloNumVar &delta_l, bool enabled);
	void appendSizeObjective();
	bool addCoordinateConstraints(const vector<double> &coord);
	bool createModel();
	void saveSolution();

public:
	SolverConfiguration CPLEXOptions;

public:
	static const double SlackMax = 5e7;
	static const double CoordMax = 1e10;
	static const double DesiredDistanceSlackMax = 6;
	static const double DesiredOrderSlackMax = 1;

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
	
	
	vector<double> coord;
};

#endif
