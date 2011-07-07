#include <iostream>
#include <ctime>
#include <cstdlib>
#include "Configuration.h"
#include "DataStore.h"
#include "DataStoreReader.h"
#include "DPSolver.h"
#include "Helpers.h"

#include "IterativeSolver.h"
#include "GASolver.h"
#include "MinMax.h"

Configuration config;
DataStore store;
//DPSolver solver;
//MIQPSolver solver;
GASolver solver;

bool readStore(const string &fileName, DataStore &store)
{
	DataStoreReader reader;
	bool result = reader.Open(fileName) && reader.Read(store);
	reader.Close();
	return result;
}

void printSolution(const Solver &solver)
{
	int n = solver.ContigCount;
	for (int i = 0; i < n; i++)
	{
		if (!solver.U[i])
			printf("%i : unused\n", i);
		else
			printf("%i : x = %3lf, u = %i, t = %i\n", i, solver.X[i], solver.U[i], solver.T[i]);
	}
	printf("Fitness %lf\n", solver.GetObjective());
}

bool getOrientation(const string &name)
{
	return name[0] == '+';
}

void testFormulation()
{
	vector<bool> uBest(store.ContigCount, 1), tBest(store.ContigCount);
	for (int i = 0; i < store.ContigCount; i++)
		tBest[i] = !getOrientation(store[i].Sequence.Name());
	IterativeSolver solver(uBest, tBest, store.ContigCount);
	solver.Options = config.Options;
	if (!solver.Formulate(store))
	{
		cout << "Unable to formulate!" << endl;
		return;
	}
	if (!solver.Solve())
	{
		cout << "Unable to solve!" << endl;
		cout << "Cplex : " << solver.GetCplexStatus() << endl;
		return;
	}
	printf("Objective: %.5lf\n", solver.GetObjective());
	printf("Original:  %.5lf\n", solver.GetHeuristicObjective());
	printf("Disabled: %i\n", solver.Disabled);
}

void getFixed(const vector<bool> &sT, const vector<bool> &bT)
{
	vector<bool> u(store.ContigCount, true);
	FixedMIQPSolver s1(u, sT, u.size()), s2(u, bT, u.size());
	s1.Options = s2.Options = config.Options;
	s1.Formulate(store), s2.Formulate(store);
	s1.Solve(), s2.Solve();
	printf("Found fixed: %.6lf\n", s1.GetObjective());
	printf("Suggested fixed: %.6lf\n", s2.GetObjective());
}

void getIterative(const vector<bool> &sT, const vector<bool> &bT)
{
	vector<bool> u(store.ContigCount, true);
	IterativeSolver i1(u, sT, u.size()), i2(u, bT, u.size());
	i1.Options = i2.Options = config.Options;
	i1.Formulate(store), i2.Formulate(store);
	i1.Solve(), i2.Solve();
	printf("Found iterative: %.6lf\n", i1.GetObjective());
	printf("Suggested iterative: %.6lf\n", i2.GetObjective());
}

void scaffolAnalyze(vector<bool> t1, vector<bool> t2)
{
	vector< vector<int> > cc;
	DPGraph graph(store);
	graph.FindConnectedComponents(cc);
	printf("Input scaffolds: %i\n", (int)cc.size());
}

void testGA(const GASolver &solver, const string &fileName)
{
	int forwardMismatch = 0, reverseMismatch = 0;
	vector<bool> tBest(store.ContigCount);
	for (int i = 0; i < store.ContigCount; i++)
	{
		tBest[i] = !getOrientation(store[i].Sequence.Name());
		if (tBest[i] != solver.T[i])
			forwardMismatch++;
		if (tBest[i] != !solver.T[i])
			reverseMismatch++;
	}
	int mismatch = min(forwardMismatch, reverseMismatch);
	GAIndividual ind(tBest, solver.matrix);
	scaffolAnalyze(solver.T, tBest);
	printf("Found: %.6lf\n", solver.GetObjective());
	printf("Suggested: %.6lf\n", ind.GetObjective());
	printf("Mismatch: %i (%i %i)\n", mismatch, forwardMismatch, reverseMismatch);
	getFixed(solver.T, tBest);
	getIterative(solver.T, tBest);
}

int countLinks(const DataStore &store)
{
	int cnt = 0;
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
		cnt++;
	return cnt;
}

void writeLog(const string &fileName, int success, double time, int linksBefore, int linksAfter)
{
	FILE *out = fopen(fileName.c_str(), "a");
	fprintf(out, "%i %.6lf %i %i\n", success, time, linksBefore, linksAfter);
	fclose(out);
}

int main(int argc, char *argv[])
{
	Helpers::ElapsedTimers.AddTimer();
	srand((unsigned int)time(NULL));
	if (config.ProcessCommandLine(argc, argv))
	{
		solver.Options = config.Options;
		if (!readStore(config.InputFileName, store))
		{
			cerr << "[-] Unable to read optimization problem (" << config.InputFileName << ")." << endl;
			return -1;
		}
		cerr << "[+] Read optimization problem (" << config.InputFileName << ")." << endl;
		if (config.PrintMatrix)
		{
			cerr << "[i] Original matrix:" << endl;
			Helpers::PrintDataStore(store);
		}
		if (config.Sort && !config.Bundle)
		{
			store.Sort();
			cerr << "[i] Sorted contig links." << endl;
		}
		else if (config.Bundle)
		{
			store.Bundle(config.Sort, config.BundlePerGroup, config.BundleAmbiguous, config.BundleDistance);
			cerr << "[i] Bundled contig links." << endl;
		}
		if ((config.Sort || config.Bundle) && config.PrintMatrix)
		{
			cerr << "[i] Optimized matrix:" << endl;
			Helpers::PrintDataStore(store);
		}
		if (!solver.Formulate(store))
		{
			cerr << "[-] Unable to formulate the optimization problem." << endl;
			return -2;
		}
		cerr << "[+] Formulated the optimization problem." << endl;
		if (!solver.Solve())
		{
			cerr << "[-] Unable to solve the optimization problem." << endl;
			return -3;
		}
		printSolution(solver);
		testGA(solver, config.OutputFileName);
		return 0;
	}
	cerr << config.LastError;
	return -1;
}
