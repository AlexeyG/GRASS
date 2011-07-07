#include "MIQPSolver.h"
#include <cstdlib>
#include "Globals.h"
#include "Helpers.h"
#include <cmath>

ILOMIPINFOCALLBACK3(BestObjectiveReachedCallback,
                    IloNum, relGap,
					IloBool, aborted,
					IloNum, bestObjective)
{
	if (!aborted  &&  hasIncumbent())
	{
		IloNum objective = getIncumbentObjValue();
		if (fabs(objective - bestObjective) <= fabs(bestObjective) * relGap)
			abort();
	}
};

MIQPSolver::MIQPSolver()
	: model(environment), x(environment), u(environment), t(environment), xi_f(environment), xi_r(environment), delta_f(environment), delta_r(environment), xi_f_ut(environment), xi_f_tt(environment), xi_r_tu(environment), xi_r_tt(environment), xi_f_uu(environment), xi_f_tu(environment), delta_f_ut(environment), delta_f_tt(environment), delta_r_tu(environment), delta_r_tt(environment), delta_f_uu(environment), delta_f_tu(environment), constraints(environment), g(environment), h(environment), p(environment), s(environment)
{
	bestObjective = 1;
	ContigCount = 0;
	status = Clean;
}

MIQPSolver::~MIQPSolver()
{
	environment.end();
}

bool MIQPSolver::Formulate(const DataStore &store)
{
	if (status != Clean)
		return false;
	if (!addContigs(store))
		return false;
	if (!addLinks(store))
		return false;
	if (!appendSizeObjective())
		return false;
	if (!createModel())
		return false;
	status = Formulated;
	return true;
}

bool MIQPSolver::Solve()
{
	if (status < Formulated)
		return false;
	if (ContigCount == 1)
	{
		U[0] = true;
		T[0] = false;
		X[0] = 0;
	}
	else
	{
		try
		{
			if (Options.UseOpportunisticSearch)
				cplex.setParam(cplex.ParallelMode, cplex.Opportunistic);
			cplex.setParam(cplex.Threads, Options.Threads);
			if (Options.SuppressOutput)
			{
				cplex.setOut(environment.getNullStream());
				//cplex.setError(environment.getNullStream());
				//cplex.setWarning(environment.getNullStream());
			}
			if (Options.TimeLimit > 0)
				cplex.setParam(cplex.TiLim, Options.TimeLimit);
			//cplex.setParam(cplex.ClockType, 1); // CPU Time
			if (Options.UseObjectiveHeuristic)
				cplex.use(BestObjectiveReachedCallback(environment, cplex.getParam(cplex.EpGap), false, bestObjective));
			//cplex.setParam(cplex.NumericalEmphasis, 1);
			if (!cplex.solve())
			{
				cout << cplex.getStatus() << endl; // REMOVE ME
				return false;
			}
			saveSolution();
		}
		catch (...)
		{
			status = Fail;
			return false;
		}
	}
	status = Success;
	return true;
}

SolverStatus MIQPSolver::GetStatus() const
{
	return status;
}

double MIQPSolver::GetObjective() const
{
	if (status == Success)
		return cplex.getObjValue();
	return -Helpers::Inf;
}

IloAlgorithm::Status MIQPSolver::GetCplexStatus() const
{
	return cplex.getStatus();
}

bool MIQPSolver::addContigs(const DataStore &store)
{
	ContigCount = store.ContigCount;
	len.resize(ContigCount);
	U.resize(ContigCount);
	T.resize(ContigCount);
	X.resize(ContigCount);
	optimized.resize(ContigCount);
	for (int i = 0; i < ContigCount; i++)
		if (!addContig(store[i]))
			return false;
	return true;
}

bool MIQPSolver::addContig(const Contig &contig)
{
	int id = contig.GetID();
	try
	{
		len[id] = contig.Sequence.Nucleotides.length();
		optimized[id] = false;
		x.add(IloNumVar(environment, 0, CoordMax));
		u.add(IloBoolVar(environment));
		t.add(IloBoolVar(environment));
		constraints.add(t[id] - u[id] <= 0);
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool MIQPSolver::addLinks(const DataStore &store)
{
	int num = 0;
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
	{
		if (!addLink(num, it->first.first, it->first.second, it->second))
			return false;
		++num;
	}
	return true;
}

bool MIQPSolver::addLink(int num, int a, int b, const ContigLink &link)
{
	// See paper for formulae
	bool e = link.EqualOrientation;
	bool r = link.ForwardOrder;
	double mu = link.Mean;
	double sigma = link.Std;
	double w = link.Weight;
	bestObjective += w;
	optimized[a] = optimized[b] = true;

	if (!addDistanceConstraint(a, b, e, r, sigma, mu))
		return false;
	if (!addOrderConstraint(a, b, e, r))
		return false;
	if (!appendOrientationObjective(a, b, e, w))
		return false;
	if (!appendDistanceObjective(a, b, e, w, xi_f[num], xi_r[num]))
		return false;
	if (!appendOrderObjective(a, b, e, w, delta_f[num], delta_r[num]))
		return false;
	return true;
}

bool MIQPSolver::addDistanceConstraint(int a, int b, bool e, bool r, double sigma, double mu)
{
	try
	{
		IloNumVar xi_f_l(environment, 0, SlackMax), xi_r_l(environment, 0, SlackMax);
		xi_f.add(xi_f_l);
		xi_r.add(xi_r_l);
		if (!e && !r)
		{
			constraints.add(x[a] - x[b] - sigma * xi_f_l <=   DistanceStdDev * sigma + mu - len[a] - len[b]);
			constraints.add(x[a] - x[b] + sigma * xi_f_l >= - DistanceStdDev * sigma + mu - len[a] - len[b]);
			constraints.add(x[b] - x[a] - sigma * xi_r_l <=   DistanceStdDev * sigma + mu - len[a] - len[b]);
			constraints.add(x[b] - x[a] + sigma * xi_r_l >= - DistanceStdDev * sigma + mu - len[a] - len[b]);
		}
		else if (!e && r)
		{
			constraints.add(x[b] - x[a] - sigma * xi_f_l <=   DistanceStdDev * sigma + mu);
			constraints.add(x[b] - x[a] + sigma * xi_f_l >= - DistanceStdDev * sigma + mu);
			constraints.add(x[a] - x[b] - sigma * xi_r_l <=   DistanceStdDev * sigma + mu);
			constraints.add(x[a] - x[b] + sigma * xi_r_l >= - DistanceStdDev * sigma + mu);
		}
		else if (e && !r)
		{
			constraints.add(x[a] - x[b] - sigma * xi_f_l <=   DistanceStdDev * sigma + mu - len[a]);
			constraints.add(x[a] - x[b] + sigma * xi_f_l >= - DistanceStdDev * sigma + mu - len[a]);
			constraints.add(x[b] - x[a] - sigma * xi_r_l <=   DistanceStdDev * sigma + mu - len[a]);
			constraints.add(x[b] - x[a] + sigma * xi_r_l >= - DistanceStdDev * sigma + mu - len[a]);
		}
		else if (e && r)
		{
			constraints.add(x[b] - x[a] - sigma * xi_f_l <=   DistanceStdDev * sigma + mu - len[b]);
			constraints.add(x[b] - x[a] + sigma * xi_f_l >= - DistanceStdDev * sigma + mu - len[b]);
			constraints.add(x[a] - x[b] - sigma * xi_r_l <=   DistanceStdDev * sigma + mu - len[b]);
			constraints.add(x[a] - x[b] + sigma * xi_r_l >= - DistanceStdDev * sigma + mu - len[b]);
		}
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool MIQPSolver::addOrderConstraint(int a, int b, bool e, bool r)
{
	try
	{
		IloNumVar delta_f_l(environment, 0, SlackMax), delta_r_l(environment, 0, SlackMax);
		delta_f.add(delta_f_l);
		delta_r.add(delta_r_l);
		if (!e && !r)
		{
			constraints.add(x[a] - x[b] + delta_f_l >= -len[b]);
			constraints.add(x[b] - x[a] + delta_r_l >= -len[a]);
		}
		else if (!e && r)
		{
			constraints.add(x[b] - x[a] + delta_f_l >= len[b]);
			constraints.add(x[a] - x[b] + delta_r_l >= len[a]);
		}
		else if (e && !r)
		{
			constraints.add(x[a] - x[b] + delta_f_l >= 0);
			constraints.add(x[b] - x[a] + delta_r_l >= len[b] - len[a]);
		}
		else if (e && r)
		{
			constraints.add(x[b] - x[a] + delta_f_l >= 0);
			constraints.add(x[a] - x[b] + delta_r_l >= len[a] - len[b]);
		}
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool MIQPSolver::appendOrientationObjective(int a, int b, bool e, double w)
{
	try
	{
		if (!e)
			g = g - w * (u[a] - 2 * t[a]) * (u[b] - 2 * t[b]);
		else
			g = g + w * (u[a] - 2 * t[a]) * (u[b] - 2 * t[b]);
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool MIQPSolver::appendDistanceObjective(int a, int b, bool e, double w, const IloNumVar &xi_f_l, const IloNumVar &xi_r_l)
{
	try
	{
		/*IloNumVar xi_f_ut_l(environment, 0, SlackMax), xi_f_tt_l(environment, 0, SlackMax), xi_r_tu_l(environment, 0, SlackMax), xi_r_tt_l(environment, 0, SlackMax), xi_f_uu_l(environment, 0, SlackMax), xi_f_tu_l(environment, 0, SlackMax);
		xi_f_uu.add(xi_f_uu_l);
		constraints.add(xi_f_uu_l - xi_f_l <= 0);
		constraints.add(xi_f_uu_l - SlackMax * u[a] <= 0);
		constraints.add(xi_f_uu_l - SlackMax * u[b] <= 0);
		constraints.add(xi_f_l - xi_f_uu_l + SlackMax * (u[a] + u[b])  <= 2 * SlackMax);
		xi_f_ut.add(xi_f_ut_l);
		constraints.add(xi_f_ut_l - xi_f_l <= 0);
		constraints.add(xi_f_ut_l - SlackMax * u[a] <= 0);
		constraints.add(xi_f_ut_l - SlackMax * t[b] <= 0);
		constraints.add(xi_f_l - xi_f_ut_l + SlackMax * (u[a] + t[b])  <= 2 * SlackMax);
		xi_f_tu.add(xi_f_tu_l);
		constraints.add(xi_f_tu_l - xi_f_l <= 0);
		constraints.add(xi_f_tu_l - SlackMax * t[a] <= 0);
		constraints.add(xi_f_tu_l - SlackMax * u[b] <= 0);
		constraints.add(xi_f_l - xi_f_tu_l + SlackMax * (t[a] + u[b])  <= 2 * SlackMax);
		xi_f_tt.add(xi_f_tt_l);
		constraints.add(xi_f_tt_l - xi_f_l <= 0);
		constraints.add(xi_f_tt_l - SlackMax * t[a] <= 0);
		constraints.add(xi_f_tt_l - SlackMax * t[b] <= 0);
		constraints.add(xi_f_l - xi_f_tt_l + SlackMax * (t[a] + t[b])  <= 2 * SlackMax);
		xi_r_tu.add(xi_r_tu_l);
		constraints.add(xi_r_tu_l - xi_r_l <= 0);
		constraints.add(xi_r_tu_l - SlackMax * t[a] <= 0);
		constraints.add(xi_r_tu_l - SlackMax * u[b] <= 0);
		constraints.add(xi_r_l - xi_r_tu_l + SlackMax * (t[a] + u[b])  <= 2 * SlackMax);
		xi_r_tt.add(xi_r_tt_l);
		constraints.add(xi_r_tt_l - xi_r_l <= 0);
		constraints.add(xi_r_tt_l - SlackMax * t[a] <= 0);
		constraints.add(xi_r_tt_l - SlackMax * t[b] <= 0);
		constraints.add(xi_r_l - xi_r_tt_l + SlackMax * (t[a] + t[b])  <= 2 * SlackMax);*/
		IloNumVar xi_f_ut_l(environment, 0, SlackMax), xi_f_tt_l(environment, 0, SlackMax), xi_r_tu_l(environment, 0, SlackMax);
		IloNumVar xi_r_tt_l(environment, 0, SlackMax), xi_f_uu_l(environment, 0, SlackMax), xi_f_tu_l(environment, 0, SlackMax);
		xi_f_ut.add(xi_f_ut_l);
		xi_f_tt.add(xi_f_tt_l);
		xi_r_tu.add(xi_r_tu_l);
		xi_r_tt.add(xi_r_tt_l);
		xi_f_uu.add(xi_f_uu_l);
		xi_f_tu.add(xi_f_tu_l);
		constraints.add(xi_f_ut_l - xi_f_l <= 0);
		constraints.add(xi_f_tt_l - xi_f_l <= 0);
		constraints.add(xi_r_tu_l - xi_r_l <= 0);
		constraints.add(xi_r_tt_l - xi_r_l <= 0);
		constraints.add(xi_f_uu_l - xi_f_l <= 0);
		constraints.add(xi_f_tu_l - xi_f_l <= 0);
		constraints.add(xi_f_ut_l - SlackMax * u[a] <= 0);
		constraints.add(xi_f_tt_l - SlackMax * t[a] <= 0);
		constraints.add(xi_r_tu_l - SlackMax * t[a] <= 0);
		constraints.add(xi_r_tt_l - SlackMax * t[a] <= 0);
		constraints.add(xi_f_uu_l - SlackMax * u[a] <= 0);
		constraints.add(xi_f_tu_l - SlackMax * t[a] <= 0);
		constraints.add(xi_f_ut_l - SlackMax * t[b] <= 0);
		constraints.add(xi_f_tt_l - SlackMax * t[b] <= 0);
		constraints.add(xi_r_tu_l - SlackMax * u[b] <= 0);
		constraints.add(xi_r_tt_l - SlackMax * t[b] <= 0);
		constraints.add(xi_f_uu_l - SlackMax * u[b] <= 0);
		constraints.add(xi_f_tu_l - SlackMax * u[b] <= 0);
		constraints.add(xi_f_l - xi_f_ut_l + SlackMax * (u[a] + t[b]) <= 2 * SlackMax);
		constraints.add(xi_f_l - xi_f_tt_l + SlackMax * (t[a] + t[b]) <= 2 * SlackMax);
		constraints.add(xi_r_l - xi_r_tu_l + SlackMax * (t[a] + u[b]) <= 2 * SlackMax);
		constraints.add(xi_r_l - xi_r_tt_l + SlackMax * (t[a] + t[b]) <= 2 * SlackMax);
		constraints.add(xi_f_l - xi_f_uu_l + SlackMax * (u[a] + u[b]) <= 2 * SlackMax);
		constraints.add(xi_f_l - xi_f_tu_l + SlackMax * (t[a] + u[b]) <= 2 * SlackMax);
		if (!e)
			h = h + w * (xi_f_ut_l - xi_f_tt_l + xi_r_tu_l - xi_r_tt_l);
		else
			h = h - w * (xi_f_ut_l - xi_f_tt_l - xi_r_tt_l) + w * (xi_f_uu_l - xi_f_tu_l);
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool MIQPSolver::appendOrderObjective(int a, int b, bool e, double w, const IloNumVar &delta_f_l, const IloNumVar &delta_r_l)
{
	try
	{
		IloNumVar delta_f_ut_l(environment, 0, SlackMax), delta_f_tt_l(environment, 0, SlackMax), delta_r_tu_l(environment, 0, SlackMax);
		IloNumVar delta_r_tt_l(environment, 0, SlackMax), delta_f_uu_l(environment, 0, SlackMax), delta_f_tu_l(environment, 0, SlackMax);
		delta_f_ut.add(delta_f_ut_l);
		delta_f_tt.add(delta_f_tt_l);
		delta_r_tu.add(delta_r_tu_l);
		delta_r_tt.add(delta_r_tt_l);
		delta_f_uu.add(delta_f_uu_l);
		delta_f_tu.add(delta_f_tu_l);
		constraints.add(delta_f_ut_l - delta_f_l <= 0);
		constraints.add(delta_f_tt_l - delta_f_l <= 0);
		constraints.add(delta_r_tu_l - delta_r_l <= 0);
		constraints.add(delta_r_tt_l - delta_r_l <= 0);
		constraints.add(delta_f_uu_l - delta_f_l <= 0);
		constraints.add(delta_f_tu_l - delta_f_l <= 0);
		constraints.add(delta_f_ut_l - SlackMax * u[a] <= 0);
		constraints.add(delta_f_tt_l - SlackMax * t[a] <= 0);
		constraints.add(delta_r_tu_l - SlackMax * t[a] <= 0);
		constraints.add(delta_r_tt_l - SlackMax * t[a] <= 0);
		constraints.add(delta_f_uu_l - SlackMax * u[a] <= 0);
		constraints.add(delta_f_tu_l - SlackMax * t[a] <= 0);
		constraints.add(delta_f_ut_l - SlackMax * t[b] <= 0);
		constraints.add(delta_f_tt_l - SlackMax * t[b] <= 0);
		constraints.add(delta_r_tu_l - SlackMax * u[b] <= 0);
		constraints.add(delta_r_tt_l - SlackMax * t[b] <= 0);
		constraints.add(delta_f_uu_l - SlackMax * u[b] <= 0);
		constraints.add(delta_f_tu_l - SlackMax * u[b] <= 0);
		constraints.add(delta_f_l - delta_f_ut_l + SlackMax * (u[a] + t[b]) <= 2 * SlackMax);
		constraints.add(delta_f_l - delta_f_tt_l + SlackMax * (t[a] + t[b]) <= 2 * SlackMax);
		constraints.add(delta_r_l - delta_r_tu_l + SlackMax * (t[a] + u[b]) <= 2 * SlackMax);
		constraints.add(delta_r_l - delta_r_tt_l + SlackMax * (t[a] + t[b]) <= 2 * SlackMax);
		constraints.add(delta_f_l - delta_f_uu_l + SlackMax * (u[a] + u[b]) <= 2 * SlackMax);
		constraints.add(delta_f_l - delta_f_tu_l + SlackMax * (t[a] + u[b]) <= 2 * SlackMax);
		if (!e)
			p = p + w * (delta_f_ut_l - delta_f_tt_l + delta_r_tu_l - delta_r_tt_l);
		else
			p = p - w * (delta_f_ut_l - delta_f_tt_l - delta_r_tt_l) + w * (delta_f_uu_l - delta_f_tu_l);
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool MIQPSolver::appendSizeObjective()
{
	try 
	{
		for (int i = 0; i < ContigCount; i++)
			s = s + u[i];
		s = s / ContigCount;
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool MIQPSolver::createModel()
{
	try
	{
		model.add(IloMaximize(environment, g - h - p + s));
		model.add(constraints);
		cplex = IloCplex(model);
	}
	catch (...)
	{
		return false;
	}
	return true;
}

void MIQPSolver::saveSolution()
{
	double minX = Helpers::Inf;
	for (int i = 0; i < ContigCount; i++)
	{
		U[i] = cplex.getValue(u[i]) == 1;
		T[i] = cplex.getValue(t[i]) == 1;
		X[i] = (U[i] && optimized[i] ? cplex.getValue(x[i]) : 0);
		if (U[i])
			minX = min(minX, (T[i] == 1 ? X[i] - (len[i] == 0 ? 0 : len[i] - 1) : X[i]));
	}
	for (int i = 0; i < ContigCount; i++)
		if (U[i]) X[i] -= minX;
}
