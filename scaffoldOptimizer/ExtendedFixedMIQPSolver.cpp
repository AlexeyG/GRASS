#include "ExtendedFixedMIQPSolver.h"
#include "Globals.h"
#include "Helpers.h"

ExtendedFixedMIQPSolver::ExtendedFixedMIQPSolver(const vector<bool> &u, const vector<bool> &t, int length)
	: model(environment), x(environment), xi(environment), delta(environment), constraints(environment), h(environment), p(environment)
{
	g = 0, s = 0;
	bestObjective = 1;
	status = Clean;
	ContigCount = length;
	U.insert(U.begin(), u.begin(), u.end());
	T.insert(T.begin(), t.begin(), t.end());
	X.resize(ContigCount);
	len.resize(ContigCount);
	optimized.resize(ContigCount, false);
}

ExtendedFixedMIQPSolver::~ExtendedFixedMIQPSolver()
{
	environment.end();
}

bool ExtendedFixedMIQPSolver::Formulate(const DataStore &store)
{
	return Formulate(store, vector<bool>(), vector<bool>());
}

bool ExtendedFixedMIQPSolver::Formulate(const DataStore &store, const vector<bool> &enabledDistance, const vector<bool> &enabledOrder, const vector<double> &coord)
{
	if (status != Clean)
		return false;
	if (!formulate(store, enabledDistance, enabledOrder) || !addCoordinateConstraints(coord) || !createModel())
	{
		status = Fail;
		return false;
	}
	status = Formulated;
	return true;
}

bool ExtendedFixedMIQPSolver::Formulate(const DataStore &store, const vector<bool> &enabledDistance, const vector<bool> &enabledOrder)
{
	if (status != Clean)
		return false;
	if (!formulate(store, enabledDistance, enabledOrder) || !createModel())
	{
		status = Fail;
		return false;
	}
	status = Formulated;
	return true;
}

bool ExtendedFixedMIQPSolver::Solve()
{
	if (status < Formulated)
		return false;
	try
	{
		cplex.setParam(cplex.ParallelMode, (Options.UseOpportunisticSearch ? -1 : 1));
		cplex.setParam(cplex.Threads, Options.LPThreads);
		if (Options.SuppressOutput)
		{
			cplex.setOut(environment.getNullStream());
			cplex.setError(environment.getNullStream());
			cplex.setWarning(environment.getNullStream());
		}
		if (Options.LPTimeLimit > 0)
			cplex.setParam(cplex.TiLim, Options.LPTimeLimit);
		//cplex.setParam(cplex.RootAlg, IloCplex::Primal);
		//cplex.setParam(cplex.PreInd, 0);
		//cplex.setParam(cplex.EpGap, 1e-8);
		cplex.setParam(cplex.NumericalEmphasis, 1);
		for (int att = 0; att < Options.LPAttempts; att++)
			if (cplex.solve())
				break;
		if (cplex.getStatus() != IloAlgorithm::Optimal)
			return false;
		saveSolution();
	}
	/*catch (IloException ex)
	{
		cout << "status: " << cplex.getStatus() << endl;
		ex.print(cout);
		return false;
	}*/
	catch (...)
	{
		status = Fail;
		return false;
	}
	status = Success;
	return true;
}

SolverStatus ExtendedFixedMIQPSolver::GetStatus() const
{
	return status;
}

double ExtendedFixedMIQPSolver::GetObjective() const
{
	if (status == Success)
		return cplex.getObjValue();
	return -Helpers::Inf;
}

IloAlgorithm::Status ExtendedFixedMIQPSolver::GetCplexStatus() const
{
	return cplex.getStatus();
}

double ExtendedFixedMIQPSolver::GetDistanceSlack(int i) const
{
	if (status == Success)
		return cplex.getValue(xi[i]);
	return -Helpers::Inf;
}

int ExtendedFixedMIQPSolver::GetSlackCount() const
{
	return xi.getSize();
}

double ExtendedFixedMIQPSolver::GetOrderSlack(int i) const
{
	if (status == Success)
		return cplex.getValue(delta[i]);
	return -Helpers::Inf;
}

bool ExtendedFixedMIQPSolver::formulate(const DataStore &store, const vector<bool> &enabledDistance, const vector<bool> &enabledOrder)
{
	if (store.ContigCount != ContigCount)
		return false;
	if (!addContigs(store))
		return false;
	if (!addLinks(store, enabledDistance, enabledOrder))
		return false;
	appendSizeObjective();
	return true;
}

bool ExtendedFixedMIQPSolver::addContigs(const DataStore &store)
{
	for (int i = 0; i < ContigCount; i++)
		if (!addContig(store[i]))
			return false;
	return true;
}

bool ExtendedFixedMIQPSolver::addContig(const Contig &contig)
{
	int id = contig.GetID();
	try
	{
		len[id] = contig.Sequence.Nucleotides.length();
		x.add(IloNumVar(environment, 0, CoordMax));
	}
	catch (IloException ex)
	{
		ex.print(cout);
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool ExtendedFixedMIQPSolver::addLinks(const DataStore &store, const vector<bool> &enabledDistance, const vector<bool> &enabledOrder)
{
	int num = 0;
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
		if (!addLink(it->first.first, it->first.second, it->second, num, (int)enabledDistance.size() <= num || enabledDistance[num], (int)enabledOrder.size() <= num || enabledOrder[num]))
			return false;
	return true;
}

bool ExtendedFixedMIQPSolver::addLink(int a, int b, const ContigLink &link, int &num, bool enabledDistance, bool enabledOrder)
{
	bool e = link.EqualOrientation;
	double w = link.Weight;
	bestObjective += w;

	if (!U[a] || !U[b])
		return true;

	if ((T[a] ^ T[b]) == e)
		return true;

	IloNumVar xi_l, delta_l;
	bool r = link.ForwardOrder;
	double mu = link.Mean;
	double sigma = link.Std;
	optimized[a] = optimized[b] = true;

	appendOrientationObjective(a, b, e, w);
	if (!addDistanceConstraint(a, b, e, r, sigma, mu, xi_l))
		return false;
	if (!addOrderConstraint(a, b, e, r, delta_l))
		return false;
	if (!appendDistanceObjective(a, b, e, w, xi_l, enabledDistance))
		return false;
	if (!appendOrderObjective(a, b, e, w, delta_l, enabledOrder))
		return false;
	num++;
	return true;
}

bool ExtendedFixedMIQPSolver::addDistanceConstraint(int a, int b, bool e, bool r, double sigma, double mu, IloNumVar &xi_l)
{
	try
	{
		xi_l = IloNumVar(environment, 0, SlackMax);
		xi.add(xi_l);
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

bool ExtendedFixedMIQPSolver::addOrderConstraint(int a, int b, bool e, bool r, IloNumVar &delta_l)
{
	try
	{
		delta_l = IloNumVar(environment, 0, SlackMax);
		delta.add(delta_l);
		if (!e && !r)
		{
			if (!T[a])
                            //constraints.add(x[a] - x[b] + delta_l * len[b] >= 0);
				constraints.add(x[a] - x[b] + delta_l * (len[b] + len[a]) / 2.0 >= 0);
			else
                            //constraints.add(x[b] - x[a] + delta_l * len[a] >= 0);
				constraints.add(x[b] - x[a] + delta_l * (len[a] + len[b]) / 2.0 >= 0);
		}
		else if (!e && r)
		{
			if (!T[a])
                            //constraints.add(x[b] - len[b] - x[a] - len[a] + delta_l * len[a] >= 0);
				constraints.add(x[b] - len[b] - x[a] - len[a] + delta_l * (len[a] + len[b]) / 2.0 >= 0);
			else
                            //constraints.add(x[a] - len[a] - x[b] - len[b] + delta_l * len[b] >= 0);
				constraints.add(x[a] - len[a] - x[b] - len[b] + delta_l * (len[b] + len[a]) / 2.0 >= 0);
		}
		else if (e && !r)
		{
			if (!T[a])
                            //constraints.add(x[a] - x[b] - len[b] + delta_l * len[b] >= 0);
				constraints.add(x[a] - x[b] - len[b] + delta_l * (len[b] + len[a]) / 2.0 >= 0);
			else
                            //constraints.add(x[b] - len[b] - x[a] + delta_l * len[a] >= 0);
				constraints.add(x[b] - len[b] - x[a] + delta_l * (len[a] + len[b]) / 2.0 >= 0);
		}
		else if (e && r)
		{
			if (!T[a])
                            //constraints.add(x[b] - x[a] - len[a] + delta_l * len[a] >= 0);
				constraints.add(x[b] - x[a] - len[a] + delta_l * (len[a] + len[b]) / 2.0 >= 0);
			else
                            //constraints.add(x[a] - len[a] - x[b] + delta_l * len[b] >= 0);
				constraints.add(x[a] - len[a] - x[b] + delta_l * (len[b] + len[a]) / 2.0 >= 0);
		}
	}
	catch (...)
	{
		return false;
	}
	return true;
}

void ExtendedFixedMIQPSolver::appendOrientationObjective(int a, int b, bool e, double w)
{
	g += w;
}

bool ExtendedFixedMIQPSolver::appendDistanceObjective(int a, int b, bool e, double w, const IloNumVar &xi_l, bool enabled)
{
	try
	{
		if (enabled)
			h = h + (xi_l / DesiredDistanceSlackMax) * w;
		else
			h = h + w;
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool ExtendedFixedMIQPSolver::appendOrderObjective(int a, int b, bool e, double w, const IloNumVar &delta_l, bool enabled)
{
	try
	{
		if (enabled)
			p = p + (delta_l / DesiredOrderSlackMax) * w;
		else
			p = p + w;
	}
	catch (...)
	{
		return false;
	}
	return true;
}

bool ExtendedFixedMIQPSolver::addCoordinateConstraints(const vector<double> &coord)
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

void ExtendedFixedMIQPSolver::appendSizeObjective()
{
	for (int i = 0; i < ContigCount; i++)
		s += U[i];
	s = s / (double)ContigCount;
}

bool ExtendedFixedMIQPSolver::createModel()
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

void ExtendedFixedMIQPSolver::saveSolution()
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
