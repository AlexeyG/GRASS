#ifndef _MIQPSOLVER_H
#define _MIQPSOLVER_H

#include "Solver.h"
#include "DataStore.h"
#include "SolverConfiguration.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

class BestObjectiveReachedCallback;

class MIQPSolver : public Solver
{
public:
	MIQPSolver();
	virtual ~MIQPSolver();

public:
	virtual bool Formulate(const DataStore &store);
	virtual bool Solve();
	virtual SolverStatus GetStatus() const;
	virtual double GetObjective() const;

public:
	IloAlgorithm::Status GetCplexStatus() const;

private:
	bool addContigs(const DataStore &store);
	bool addContig(const Contig &contig);
	bool addLinks(const DataStore &store);
	bool addLink(int num, int a, int b, const ContigLink &link);
	bool addDistanceConstraint(int a, int b, bool e, bool r, double sigma, double mu);
	bool addOrderConstraint(int a, int b, bool e, bool r);
	bool appendOrientationObjective(int a, int b, bool e, double w);
	bool appendDistanceObjective(int a, int b, bool e, double w, const IloNumVar &xi_f_l, const IloNumVar &xi_r_l);
	bool appendOrderObjective(int a, int b, bool e, double w, const IloNumVar &delta_f_l, const IloNumVar &delta_r_l);
	bool appendSizeObjective();
	bool createModel();
	void saveSolution();

public:
	static const double SlackMax = 1e7;
	static const double CoordMax = 1e10;
	static const double DistanceStdDev = 1.0;

private:
	IloEnv environment;
	IloModel model;
	IloNumVarArray x,u,t;
	IloNumVarArray xi_f, xi_r, delta_f, delta_r;
	IloNumVarArray xi_f_ut, xi_f_tt, xi_r_tu, xi_r_tt, xi_f_uu, xi_f_tu;
	IloNumVarArray delta_f_ut, delta_f_tt, delta_r_tu, delta_r_tt, delta_f_uu, delta_f_tu;
	IloRangeArray constraints;
	IloExpr g,h,p,s;
	IloCplex cplex;
	double bestObjective;
	vector<int> len;
	vector<bool> optimized;
};
#endif
