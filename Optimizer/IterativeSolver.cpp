#include "IterativeSolver.h"
#include "Helpers.h"

IterativeSolver::IterativeSolver(const vector<bool> &u, const vector<bool> &t, int length)
	: U(u), T(t), solver(u, t, length), extension(u, t, length)
{
	status = Clean;
	Disabled = 0;
	ContigCount = length;
	X.resize(ContigCount);
}

IterativeSolver::~IterativeSolver()
{
}

bool IterativeSolver::Formulate(const DataStore &store)
{
	if (status != Clean)
		return false;
	if (!solver.Formulate(store))
		return false;
	this->store = store;
	solver.Options = extension.Options = Options;
	status = Formulated;
	return true;
}

bool IterativeSolver::Solve()
{
	if (status != Formulated)
		return false;
	if (!solver.Solve())
	{
		status = Fail;
		return false;
	}
	int size = solver.GetSlackCount();
	vector<bool> distance(size, true), order(size, true);
	for (int i = 0; i < size; i++)
	{
		if (solver.GetDistanceSlack(i) > FixedMIQPSolver::DesiredDistanceSlackMax + Helpers::Eps)
			distance[i] = false, Disabled++;
		if (solver.GetOrderSlack(i) > FixedMIQPSolver::DesiredOrderSlackMax + Helpers::Eps)
			order[i] = false, Disabled++;
	}
	if (!extension.Formulate(store, distance, order) || !extension.Solve())
	{
		status = Fail;
		return false;
	}
	X = extension.X;
	status = Success;
	return true;
}

SolverStatus IterativeSolver::GetStatus() const
{
	return status;
}

double IterativeSolver::GetObjective() const
{
	return extension.GetObjective();
}

IloAlgorithm::Status IterativeSolver::GetCplexStatus() const
{
	return extension.GetCplexStatus();
}

double IterativeSolver::GetDistanceSlack(int i) const
{
	return extension.GetDistanceSlack(i);
}

double IterativeSolver::GetOrderSlack(int i) const
{
	return extension.GetOrderSlack(i);
}

int IterativeSolver::GetSlackCount() const
{
	return extension.GetSlackCount();
}

double IterativeSolver::GetHeuristicObjective() const
{
	return solver.GetObjective();
}
