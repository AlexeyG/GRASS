#include "DPSolver.h"
#include "MIQPSolver.h"
#include "Helpers.h"
#include "MinMax.h"

DPSolver::DPSolver()
	: nComponents(0)
{
	ContigCount = 0;
	status = Clean;
	objectiveValue = 0;
}

DPSolver::~DPSolver()
{
}

bool DPSolver::Formulate(const DataStore &store)
{
	if (status != Clean)
		return false;
	ContigCount = store.ContigCount;
	U.resize(ContigCount);
	T.resize(ContigCount);
	X.resize(ContigCount);
	this->store = store;
	graph = DPGraph(store);
	nComponents = graph.FindConnectedComponents(connectedComponents);
	status = Formulated;
	return true;
}

bool DPSolver::Solve()
{
	bool result;
	if (status < Formulated)
		return false;
	if (nComponents == 0)
		result = true;
	else if (nComponents > 1)
		result = processComponents();
	else
		result = processSingleProblem();
	status = (result ? Success : Fail);
	return result;
}

SolverStatus DPSolver::GetStatus() const
{
	return status;
}

double DPSolver::GetObjective() const
{
	if (status == Success)
		return objectiveValue;
	return -Helpers::Inf;
}

bool DPSolver::processSingleProblem()
{
	MIQPSolver solver;
	solver.Options = Options;
	int nContigs = store.ContigCount;
	bool form = solver.Formulate(store);
	bool solve = form && solver.Solve();
	if (form && solve)
	{
		for (int i = 0; i < nContigs; i++)
			U[i] = solver.U[i], T[i] = solver.T[i], X[i] = solver.X[i];
		objectiveValue = solver.GetObjective();
		return true;
	}
	return false;
}

bool DPSolver::processComponents()
{
	bool result = true;
	vector<double> minX(nComponents);
	vector<double> maxX(nComponents);
	for (int i = 0; i < nComponents; i++)
	{
		DataStore compStore;
		MIQPSolver solver;
		solver.Options = Options;
		compStore.Extract(connectedComponents[i], compStore);
		if (!solver.Formulate(compStore) || !solver.Solve())
		{
			result = false;
			break;
		}
		objectiveValue += solver.GetObjective();
		minX[i] =   Helpers::Inf;
		maxX[i] = - Helpers::Inf;
		int nContigsComponent = connectedComponents[i].size();
		for (int j = 0; j < nContigsComponent; j++)
		{
			U[connectedComponents[i][j]] = solver.U[j];
			T[connectedComponents[i][j]] = solver.T[j];
			X[connectedComponents[i][j]] = solver.X[j];
			if (solver.U[j])
			{
				minX[i] = min(minX[i], (T[i] == 1 ? solver.X[j] - compStore[j].Sequence.Nucleotides.length() + 1 : solver.X[j]));
				maxX[i] = max(maxX[i], (T[i] == 0 ? solver.X[j] + compStore[j].Sequence.Nucleotides.length() - 1 : solver.X[j]));
			}
		}
	}
	for (int i = 0; i < nComponents; i++)
	{
		int shift = minX[i];
		int offset = (i > 0 ? maxX[i - 1] + ScaffoldSeprator : 0);
		int nContigsComponent = connectedComponents[i].size();
		for (int j = 0; j < nContigsComponent; j++)
			if (U[connectedComponents[i][j]])
				X[connectedComponents[i][j]] -= shift - offset;
	}
	return result;
}
