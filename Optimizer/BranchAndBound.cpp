#include "BranchAndBound.h"
#include "Globals.h"
#include "Helpers.h"
#include "ExtendedFixedMIQPSolver.h"
#include "FixedMIQPSolver.h"

BranchAndBound::BranchAndBound(const vector<bool> &u, const vector<bool> &t, int length)
	: model(environment), x(environment), xi(environment), delta(environment), alpha(environment), beta(environment), xi_alpha(environment), delta_beta(environment), constraints(environment), h(environment), p(environment)
{
	g = 0, s = 0;
	status = Clean;
	ContigCount = length;
	U.insert(U.begin(), u.begin(), u.end());
	T.insert(T.begin(), t.begin(), t.end());
	X.resize(ContigCount);
	len.resize(ContigCount);
	optimized.resize(ContigCount, false);
}

BranchAndBound::~BranchAndBound()
{
	environment.end();
}

bool BranchAndBound::Formulate(const DataStore &store, const vector<double> &coord)
{
	if (status != Clean)
		return false;
	if (!formulate(store) || !addCoordinateConstraints(coord) || !createModel() || !assignPriorities(store))
	{
		status = Fail;
		return false;
	}
	status = Formulated;
	return true;
}

bool BranchAndBound::Formulate(const DataStore &store)
{
	if (status != Clean)
		return false;
	if (!formulate(store) || !createModel() || !assignPriorities(store))
	{
		status = Fail;
		return false;
	}
	status = Formulated;
	return true;
}

bool BranchAndBound::Solve()
{
	if (status < Formulated)
		return false;
	try
	{
		cplex.setParam(cplex.ParallelMode, Options.UseOpportunisticSearch);
		cplex.setParam(cplex.Threads, Options.Threads);
		if (Options.SuppressOutput)
		{
			cplex.setOut(environment.getNullStream());
			cplex.setError(environment.getNullStream());
			cplex.setWarning(environment.getNullStream());
		}
		//cplex.setParam(cplex.ClockType, 1); // CPU Time
		cplex.setParam(cplex.NumericalEmphasis, 1);
		cplex.setParam(cplex.MIPOrdInd, 1);
		//cplex.exportModel("model.mps");
		//cplex.writeMIPStart("start.mst");
		if (Options.TimeLimit > 0)
			cplex.setParam(cplex.TiLim, Options.TimeLimit);
		if (!cplex.solve())
			return false;
		saveSolution();
	}
	catch (IloException &ex)
	{
		ex.print(cout);
		return false;
	}
	catch (...)
	{
		status = Fail;
		return false;
	}
	status = Success;
	return true;
}

SolverStatus BranchAndBound::GetStatus() const
{
	return status;
}

double BranchAndBound::GetObjective() const
{
	if (status == Success)
		return cplex.getObjValue();
	return -Helpers::Inf;
}

IloAlgorithm::Status BranchAndBound::GetCplexStatus() const
{
	return cplex.getStatus();
}

double BranchAndBound::GetDistanceSlack(int i) const
{
	if (status == Success)
		return cplex.getValue(xi[i]);
	return -Helpers::Inf;
}

double BranchAndBound::GetOrderSlack(int i) const
{
	if (status == Success)
		return cplex.getValue(delta[i]);
	return -Helpers::Inf;
}

int BranchAndBound::GetSlackCount() const
{
	return xi.getSize();
}

bool BranchAndBound::formulate(const DataStore &store)
{
	if (store.ContigCount != ContigCount)
		return false;
	if (!addContigs(store))
		return false;
	if (!addLinks(store))
		return false;
	appendSizeObjective();
	return true;
}

bool BranchAndBound::addContigs(const DataStore &store)
{
	for (int i = 0; i < ContigCount; i++)
		if (!addContig(store[i]))
			return false;
	return true;
}

bool BranchAndBound::addContig(const Contig &contig)
{
	int id = contig.GetID();
	try
	{
		len[id] = contig.Sequence.Nucleotides.length();
		x.add(IloNumVar(environment, 0, CoordMax));
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool BranchAndBound::addLinks(const DataStore &store)
{
	int num = 0;
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
		if (!addLink(it->first.first, it->first.second, it->second, num))
			return false;
	return true;
}

bool BranchAndBound::addLink(int a, int b, const ContigLink &link, int &num)
{
	bool e = link.EqualOrientation;
	double w = link.Weight;

	if (!U[a] || !U[b])
		return true;

	if ((T[a] ^ T[b]) == e)
		return true;

	bool r = link.ForwardOrder;
	double mu = link.Mean;
	double sigma = link.Std;
	optimized[a] = optimized[b] = true;

	appendOrientationObjective(a, b, e, w);
	if (!addDistanceConstraint(a, b, e, r, sigma, mu))
		return false;
	if (!addOrderConstraint(a, b, e, r))
		return false;
	if (!appendDistanceObjective(a, b, e, w, num))
		return false;
	if (!appendOrderObjective(a, b, e, w, num))
		return false;
	num++;
	return true;
}

bool BranchAndBound::addDistanceConstraint(int a, int b, bool e, bool r, double sigma, double mu)
{
	try
	{
		IloNumVar xi_l(environment, 0, SlackMax);
		IloBoolVar alpha_l(environment);
		IloNumVar xi_alpha_l(environment, 0, SlackMax);
		xi.add(xi_l);
		alpha.add(alpha_l);
		xi_alpha.add(xi_alpha_l);

		constraints.add(xi_alpha_l - alpha_l * SlackMax <= 0);
		constraints.add(xi_alpha_l - xi_l <= 0);
		constraints.add(xi_l - xi_alpha_l + SlackMax * alpha_l <= SlackMax);
		if (!e && !r)
		{
			if (!T[a])
			{
				constraints.add(x[a] - x[b] - sigma * xi_l <=   sigma + mu - len[a] - len[b]);
				constraints.add(x[a] - x[b] + sigma * xi_l >= - sigma + mu - len[a] - len[b]);
			}
			else
			{
				constraints.add(x[b] - x[a] - sigma * xi_l <=   sigma + mu - len[a] - len[b]);
				constraints.add(x[b] - x[a] + sigma * xi_l >= - sigma + mu - len[a] - len[b]);
			}
		}
		else if (!e && r)
		{
			if (!T[a])
			{
				constraints.add(x[b] - x[a] - sigma * xi_l <=   sigma + mu);
				constraints.add(x[b] - x[a] + sigma * xi_l >= - sigma + mu);
			}
			else
			{
				constraints.add(x[a] - x[b] - sigma * xi_l <=   sigma + mu);
				constraints.add(x[a] - x[b] + sigma * xi_l >= - sigma + mu);
			}
		}
		else if (e && !r)
		{
			if (!T[a])
			{
				constraints.add(x[a] - x[b] - sigma * xi_l <=   sigma + mu - len[a]);
				constraints.add(x[a] - x[b] + sigma * xi_l >= - sigma + mu - len[a]);
			}
			else
			{
				constraints.add(x[b] - x[a] - sigma * xi_l <=   sigma + mu - len[a]);
				constraints.add(x[b] - x[a] + sigma * xi_l >= - sigma + mu - len[a]);
			}
		}
		else if (e && r)
		{
			if (!T[a])
			{
				constraints.add(x[b] - x[a] - sigma * xi_l <=   sigma + mu - len[b]);
				constraints.add(x[b] - x[a] + sigma * xi_l >= - sigma + mu - len[b]);
			}
			else
			{
				constraints.add(x[a] - x[b] - sigma * xi_l <=   sigma + mu - len[b]);
				constraints.add(x[a] - x[b] + sigma * xi_l >= - sigma + mu - len[b]);
			}
		}
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool BranchAndBound::addOrderConstraint(int a, int b, bool e, bool r)
{
	try
	{
		IloNumVar delta_l(environment, 0, SlackMax);
		IloBoolVar beta_l(environment);
		IloNumVar delta_beta_l(environment, 0, SlackMax);
		delta.add(delta_l);
		beta.add(beta_l);
		delta_beta.add(delta_beta_l);

		constraints.add(delta_beta_l - beta_l * SlackMax <= 0);
		constraints.add(delta_beta_l - delta_l <= 0);
		constraints.add(delta_l - delta_beta_l + SlackMax * beta_l <= SlackMax);
		if (!e && !r)
		{
			if (!T[a])
				//constraints.add(x[a] - x[b] + delta_l * len[b] >= -len[b]);
				constraints.add(x[a] - x[b] + delta_l * len[b] >= 0);
			else
				//constraints.add(x[b] - x[a] + delta_l * len[a] >= -len[a]);
				constraints.add(x[b] - x[a] + delta_l * len[a] >= 0);
		}
		else if (!e && r)
		{
			if (!T[a])
				//constraints.add(x[b] - x[a] + delta_l * len[a] >= len[b]);
				constraints.add(x[b] - len[b] - x[a] - len[a] + delta_l * len[a] >= 0);
			else
				//constraints.add(x[a] - x[b] + delta_l * len[b] >= len[a]);
				constraints.add(x[a] - len[a] - x[b] - len[b] + delta_l * len[b] >= 0);
		}
		else if (e && !r)
		{
			if (!T[a])
				//constraints.add(x[a] - x[b] + delta_l * len[b] >= 0);
				constraints.add(x[a] - x[b] - len[b] + delta_l * len[b] >= 0);
			else
				//constraints.add(x[b] - x[a] + delta_l * len[a] >= len[b] - len[a]);
				constraints.add(x[b] - len[b] - x[a] + delta_l * len[a] >= 0);
		}
		else if (e && r)
		{
			if (!T[a])
				//constraints.add(x[b] - x[a] + delta_l * len[a] >= 0);
				constraints.add(x[b] - x[a] - len[a] + delta_l * len[a] >= 0);
			else
				//constraints.add(x[a] - x[b] + delta_l * len[b] >= len[a] - len[b]);
				constraints.add(x[a] - len[a] - x[b] + delta_l * len[b] >= 0);
		}
	}
	catch (...)
	{
		return false;
	}
	return true;
}

void BranchAndBound::appendOrientationObjective(int a, int b, bool e, double w)
{
	g += w;
}

bool BranchAndBound::appendDistanceObjective(int a, int b, bool e, double w, int num)
{
	try
	{
		h = h + (1 - alpha[num] + xi_alpha[num] / ExtendedFixedMIQPSolver::DesiredDistanceSlackMax) * w;
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool BranchAndBound::appendOrderObjective(int a, int b, bool e, double w, int num)
{
	try
	{
		p = p + (1 - beta[num] + delta_beta[num] / ExtendedFixedMIQPSolver::DesiredOrderSlackMax) * w;
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool BranchAndBound::addCoordinateConstraints(const vector<double> &coord)
{
	if ((int)coord.size() != ContigCount)
		return false;
	try
	{
		for (int i = 0; i < ContigCount; i++)
			constraints.add(x[i] == coord[i]);
	}
	catch (...)
	{
		return false;
	}
	return true;
}

void BranchAndBound::appendSizeObjective()
{
	for (int i = 0; i < ContigCount; i++)
		s += U[i];
	s = s / (double)ContigCount;
}

bool BranchAndBound::createModel()
{
	try
	{
		model.add(IloMaximize(environment, (g + s) - 0.5 * h - 0.5 * p));
		model.add(constraints);
		cplex = IloCplex(model);
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool BranchAndBound::assignPriorities(const DataStore &store)
{
	FixedMIQPSolver solver(U, T, ContigCount);
	solver.Options = Options;
	if (!solver.Formulate(store) || !solver.Solve())
		return false;
	try
	{
		IloNumVarArray vars(environment);
		vars.add(alpha);
		IloNumArray start(environment);
		for (int i = 0; i < alpha.getSize(); i++)
		{
			double slack = solver.GetDistanceSlack(i);
			cplex.setPriority(alpha[i], slack);
			Slack.push_back(slack);
			if (slack > ExtendedFixedMIQPSolver::DesiredDistanceSlackMax)
				start.add(0), Incumbent.push_back(false);
			else
				start.add(1), Incumbent.push_back(true);
		}
		vars.add(beta);
		for (int i = 0; i < beta.getSize(); i++)
		{
			double slack = solver.GetOrderSlack(i);
			cplex.setPriority(beta[i], slack);
			Slack.push_back(slack);
			if (slack > ExtendedFixedMIQPSolver::DesiredOrderSlackMax)
				start.add(0), Incumbent.push_back(false);
			else
				start.add(1), Incumbent.push_back(true);
		}
		cplex.addMIPStart(vars, start);
	}
	catch (...)
	{
		return false;
	}
	return true;
}

void BranchAndBound::saveSolution()
{
	double minX = Helpers::Inf;
	for (int i = 0; i < ContigCount; i++)
	{
		if (!U[i])
			T[i] = false;
		X[i] = (U[i] && optimized[i] ? cplex.getValue(x[i]) : 0);
		if (U[i])
			minX = min(minX, (T[i] == 1 ? X[i] - (len[i] == 0 ? 0 : len[i] - 1) : X[i]));
	}
	for (int i = 0; i < ContigCount; i++)
		if (U[i]) X[i] -= minX;
}
