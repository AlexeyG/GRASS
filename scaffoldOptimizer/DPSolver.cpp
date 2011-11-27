#include "DPSolver.h"
#include "EMSolver.h"
#include "Helpers.h"
#include "MinMax.h"

DPSolver::DPSolver()
	: nComponents(0)
{
	ContigCount = 0;
	status = Clean;
	objectiveValue = 0;
	MaxIteration = 0;
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
	fprintf(stderr, "[i] Have %i connected components.\n", nComponents);
	if (nComponents == 0)
		result = true;
	else
		result = processComponents();
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

bool DPSolver::processComponents()
{
	bool result = true;
	vector<double> minX(nComponents);
	vector<double> maxX(nComponents);
	vector< vector<int> > backTransform(nComponents);
	vector<Scaffold> scaffolds;
	MaxIteration = 0;
	for (int i = 0; i < nComponents; i++)
	{
		int nContigsComponent = connectedComponents[i].size();
		fprintf(stderr, "    [i] Processing component %i of size %i.\n", i + 1, nContigsComponent);
		DataStore compStore;
		EMSolver *solver = new EMSolver();
		solver->Options = Options;
		store.Extract(connectedComponents[i], compStore, backTransform[i]);
		if (!solver->Formulate(compStore) || !solver->Solve())
		{
			fprintf(stderr, "        [-] Unable to solve or formulate.\n");
			result = false;
			delete solver;
			break;
		}
		else
			fprintf(stderr, "        [+] Formulated and solved subproblem.\n");
		scaffolds = ScaffoldExtractor::Extract(*solver);
		for (vector<Scaffold>::iterator it = scaffolds.begin(); it != scaffolds.end(); it++)
		{
			it->ApplyTransform(backTransform[i]);
			it->NormalizeCoordindates();
		}
		Scaffolds.insert(Scaffolds.end(), scaffolds.begin(), scaffolds.end());
		objectiveValue += solver->GetObjective();
		MaxIteration = max(MaxIteration, solver->Iteration);
		minX[i] =   Helpers::Inf;
		maxX[i] = - Helpers::Inf;
		for (int j = 0; j < nContigsComponent; j++)
		{
			int id = backTransform[i][j];
			U[id] = solver->U[j];
			T[id] = solver->T[j];
			X[id] = solver->X[j];
			if (solver->U[j])
			{
				int contigLen = compStore[j].Sequence.Nucleotides.length();
				minX[i] = min(minX[i], (solver->T[j] == 1 ? solver->X[j] - contigLen + 1 : solver->X[j]));
				maxX[i] = max(maxX[i], (solver->T[j] == 0 ? solver->X[j] + contigLen - 1 : solver->X[j]));
			}
		}
		delete solver;
	}
	fprintf(stderr, "    [i] Adjusting subscaffold positions\n");
	for (int i = 0; i < nComponents; i++)
	{
		int shift = minX[i];
		int offset = (i > 0 ? maxX[i - 1] + ScaffoldSeprator : 0); // that's a strange statement - was I planning to put all scaffolds on a single line?
		//int offset = 0;
		int nContigsComponent = connectedComponents[i].size();
		for (int j = 0; j < nContigsComponent; j++)
			if (U[backTransform[i][j]])
				X[backTransform[i][j]] -= shift - offset;
	}
	return result;
}
