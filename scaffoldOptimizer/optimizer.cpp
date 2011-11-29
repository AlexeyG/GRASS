#include <iostream>
#include <ctime>
#include <cstdlib>
#include "Configuration.h"
#include "DataStore.h"
#include "DataStoreReader.h"
#include "ReadCoverageReader.h"
#include "DPSolver.h"
#include "Helpers.h"

#include "ScaffoldExtractor.h"
#include "ScaffoldConverter.h"

#include "Writer.h"

#ifdef _INHOUSETESTS
#include "IterativeSolver.h"
#include "GASolver.h"
#include "MinMax.h"
#include "BranchAndBound.h"
#include "EMSolver.h"

#include "ScaffoldComparer.h"

#include "GraphViz.h"
#endif

using namespace std;

Configuration config;
DataStore store;
DPSolver solver;
ReadCoverage coverage;

#ifdef _INHOUSETESTS
FILE *out;

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
	return name[0] != '+';
}

int getDistance(const ScaffoldContig &a, const ScaffoldContig &b)
{
	int len = b.X - a.X;
	if (!a.T)
		len -= a.Length;
	if (b.T)
		len -= b.Length;
	return len;
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

	vector<Scaffold> i1s = ScaffoldExtractor::Extract(store, s1);
	printf("Found fixed scaffolds: %i\n", (int)i1s.size());
	for (int j = 0; j < (int)i1s.size(); j++)
	{
		printf("Scaffold %i:\n", j);
		for (int i = 0; i < i1s[j].ContigCount(); i++)
			printf("%c%i ", (i1s[j][i].T ? '-' : '+'), i1s[j][i].Id);
		printf("\n");
		for (int i = 1; i < i1s[j].ContigCount(); i++)
			printf("%i ", getDistance(i1s[j][i - 1], i1s[j][i]));
		printf("\n");
	}

	vector<Scaffold> i2s = ScaffoldExtractor::Extract(store, s2);
	printf("Suggested fixed scaffolds: %i\n", (int)i2s.size());
	for (int j = 0; j < (int)i2s.size(); j++)
	{
		printf("Scaffold %i:\n", j);
		for (int i = 0; i < i2s[j].ContigCount(); i++)
			printf("%c%i ", (i2s[j][i].T ? '-' : '+'), i2s[j][i].Id);
		printf("\n");
		for (int i = 1; i < i2s[j].ContigCount(); i++)
			printf("%i ", getDistance(i2s[j][i - 1], i2s[j][i]));
		printf("\n");
	}
	for (int i = 0; i < s2.GetSlackCount(); i++)
		printf("xi = %.3lf delta = %.3lf\n", s2.GetDistanceSlack(i), s2.GetOrderSlack(i));
}

void getIterative(const vector<bool> &sT, const vector<bool> &bT, const vector<Scaffold> &truth)
{
	vector<bool> u(store.ContigCount, true);
	IterativeSolver i1(u, sT, u.size()), i2(u, bT, u.size());
	i1.Options = i2.Options = config.Options;
	i1.Formulate(store), i2.Formulate(store);
	i1.Solve(), i2.Solve();
	printf("Found iterative: %.6lf\n", i1.GetObjective());
	printf("Suggested iterative: %.6lf\n", i2.GetObjective());

	vector<Scaffold> i1s = ScaffoldExtractor::Extract(i1), i2s = ScaffoldExtractor::Extract(i2);

	printf("Found iterative scaffolds: %i, mismatch %i\n", (int)i1s.size(), ScaffoldComparer::Compare(i1s, truth));
	for (int i = 0; i < i1s[0].ContigCount(); i++)
		printf("%c%i ", (i1s[0][i].T ? '-' : '+'), i1s[0][i].Id);
	printf("\n");
	for (int i = 1; i < i1s[0].ContigCount(); i++)
		printf("%i ", getDistance(i1s[0][i - 1], i1s[0][i]));
	printf("\n");

	printf("Suggested iterative scaffolds: %i, mismatch %i\n", (int)i2s.size(), ScaffoldComparer::Compare(i2s, truth));
	for (int i = 0; i < i2s[0].ContigCount(); i++)
		printf("%c%i ", (i2s[0][i].T ? '-' : '+'), i2s[0][i].Id);
	printf("\n");
	for (int i = 1; i < i2s[0].ContigCount(); i++)
		printf("%i ", getDistance(i2s[0][i - 1], i2s[0][i]));
	printf("\n");
	//for (int i = 0; i < i2.GetSlackCount(); i++)
	//	printf("xi = %.3lf delta = %.3lf\n", i2.GetDistanceSlack(i), i2.GetOrderSlack(i));
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

double getRealDistance(const ScaffoldContig &a, const ScaffoldContig &b)
{
	double len = 0;
	if (a.X < b.X)
	{
		len += b.X - a.X;
		if (a.T)
			len += a.Length;
		if (!b.T)
			len += b.Length;
	}
	else
	{
		len += a.X - b.X;
		if (b.T)
			len += b.Length;
		if (!a.T)
			len += a.Length;
	}
	return len;
}

ScaffoldContig getWithId(const vector<Scaffold> &s, int id)
{
	int n = s.size();
	for (int i = 0; i < n; i++)
		for (int j = 0; j < s[i].ContigCount(); j++)
			if (s[i][j].Id == id)
				return s[i][j];
	return ScaffoldContig(-1, false, 0, 0);
}

void testGround(const DPSolver &found)
{
	vector<Scaffold> truth = ScaffoldExtractor::Extract(store);
	for (int i = 0; i < (int)truth.size(); i++)
		truth[i].NormalizeCoordindates();
	vector<double> coordBest(store.ContigCount);
	vector<bool> tBest(store.ContigCount);
	for (int i = 0; i < (int)truth.size(); i++)
		for (int j = 0; j < truth[i].ContigCount(); j++)
			coordBest[truth[i][j].Id] = truth[i][j].X, tBest[truth[i][j].Id] = truth[i][j].T;
	IterativeSolver want(vector<bool>(tBest.size(), true), tBest, tBest.size());
	want.Options = config.Options;
	if (!want.Formulate(store, coordBest) || !want.Solve())
	{
		fprintf(out, "Unable to formulate or solve WANT\n");
		return;
	}
	vector<Scaffold> foundScaffold = ScaffoldExtractor::Extract(found);
	vector<Scaffold> wantScaffold = ScaffoldExtractor::Extract(want);
	
	fprintf(out, "Found scaffolds: %i\n", (int)foundScaffold.size());
	for (int j = 0; j < (int)foundScaffold.size(); j++)
	{
		for (int i = 0; i < foundScaffold[j].ContigCount(); i++)
			fprintf(out, "%c%i ", (foundScaffold[j][i].T ? '-' : '+'), foundScaffold[j][i].Id);
		fprintf(out, "\n");
		for (int i = 1; i < foundScaffold[j].ContigCount(); i++)
			fprintf(out, "%i ", getDistance(foundScaffold[j][i - 1], foundScaffold[j][i]));
		fprintf(out, "\n");
	}

	fprintf(out, "Want scaffolds: %i\n", (int)wantScaffold.size());
	for (int j = 0; j < (int)wantScaffold.size(); j++)
	{
		for (int i = 0; i < wantScaffold[j].ContigCount(); i++)
			fprintf(out, "%c%i ", (wantScaffold[j][i].T ? '-' : '+'), wantScaffold[j][i].Id);
		fprintf(out, "\n");
		for (int i = 1; i < wantScaffold[j].ContigCount(); i++)
			fprintf(out, "%i ", getDistance(wantScaffold[j][i - 1], wantScaffold[j][i]));
		fprintf(out, "\n");
	}
	fprintf(out, "Found: %.5lf Want: %.5lf\n", found.GetObjective(), want.GetObjective());
	fprintf(out, "Orientation mismatch: %i\n", ScaffoldComparer::OrientationDistance(foundScaffold, wantScaffold));
	fprintf(out, "Found distance: %i\n", ScaffoldComparer::Compare(foundScaffold, truth));
	fprintf(out, "Want distance: %i\n", ScaffoldComparer::Compare(wantScaffold, truth));
	fprintf(out, "EM iterations: %i\n", solver.MaxIteration);
	
	/*int id = 0;
    int pos = 0;
    for (DataStore::LinkMap::iterator it = store.Begin(); it != store.End(); it++, pos++)
    {
		int a = it->first.first, b = it->first.second;
		if ((tBest[a] ^ tBest[b]) != it->second.EqualOrientation)
		{
			fprintf(stderr, "%i - Linked: d(%i,%i) = %8.2f +/- %8.2f; orientation: %8s; order: %7s with weight %.5lf.\n", pos, it->second.First, it->second.Second, it->second.Mean, it->second.Std, (it->second.EqualOrientation ? "equal" : "opposite"), (it->second.ForwardOrder ? "forward" : "reverse"), it->second.Weight);
			fprintf(stderr, "Distance: %8.2f - %8.2f\n", getRealDistance(getWithId(foundScaffold, a), getWithId(foundScaffold, b)), getRealDistance(getWithId(wantScaffold, a), getWithId(wantScaffold, b)));
			printf("xi = %.3lf delta = %.3lf | xi = %.3lf delta = %.3lf\n", found.GetDistanceSlack(id), found.GetOrderSlack(id), want.GetDistanceSlack(id), want.GetOrderSlack(id));
			printf("\n");
			id++;
		}
	}*/

	/*GraphViz graphWant(wantScaffold, store, want), graphFound(foundScaffold, store, found);
	FILE *wantFile = fopen((config.OutputFileName + "_want.gv").c_str(), "w"), *foundFile = fopen((config.OutputFileName + "_found.gv").c_str(), "w");
	fprintf(wantFile, "%s\n", graphWant.GetString().c_str());
	fprintf(foundFile, "%s\n", graphFound.GetString().c_str());
	fclose(wantFile);
	fclose(foundFile);*/
}

void testGroundIterative()
{
	vector<Scaffold> truth = ScaffoldExtractor::Extract(store);
	for (int i = 0; i < (int)truth.size(); i++)
		truth[i].NormalizeCoordindates();
	vector<double> coordBest(store.ContigCount);
	vector<bool> tBest(store.ContigCount);
	for (int i = 0; i < (int)truth.size(); i++)
		for (int j = 0; j < truth[i].ContigCount(); j++)
			coordBest[truth[i][j].Id] = truth[i][j].X, tBest[truth[i][j].Id] = truth[i][j].T;
	IterativeSolver found(vector<bool>(tBest.size(), true), tBest, tBest.size());
	IterativeSolver want(vector<bool>(tBest.size(), true), tBest, tBest.size());
	found.Options = want.Options = config.Options;
	if (!found.Formulate(store) || !found.Solve())
	{
		fprintf(out, "Unable to formulate or solve FOUND\n");
		return;
	}
	if (!want.Formulate(store, coordBest) || !want.Solve())
	{
		fprintf(out, "Unable to formulate or solve WANT\n");
		return;
	}
	vector<Scaffold> foundScaffold = ScaffoldExtractor::Extract(found);
	vector<Scaffold> wantScaffold = ScaffoldExtractor::Extract(want);
	fprintf(out, "Found: %.5lf Want: %.5lf\n", found.GetObjective(), want.GetObjective());
	fprintf(out, "Found scaffolds: %i\n", (int)foundScaffold.size());
	for (int j = 0; j < (int)foundScaffold.size(); j++)
	{
		for (int i = 0; i < foundScaffold[j].ContigCount(); i++)
			fprintf(out, "%c%i ", (foundScaffold[j][i].T ? '-' : '+'), foundScaffold[j][i].Id);
		fprintf(out, "\n");
		for (int i = 1; i < foundScaffold[j].ContigCount(); i++)
			fprintf(out, "%i ", getDistance(foundScaffold[j][i - 1], foundScaffold[j][i]));
		fprintf(out, "\n");
	}

	fprintf(out, "Want scaffolds: %i\n", (int)wantScaffold.size());
	for (int j = 0; j < (int)wantScaffold.size(); j++)
	{
		for (int i = 0; i < wantScaffold[j].ContigCount(); i++)
			fprintf(out, "%c%i ", (wantScaffold[j][i].T ? '-' : '+'), wantScaffold[j][i].Id);
		fprintf(out, "\n");
		for (int i = 1; i < wantScaffold[j].ContigCount(); i++)
			fprintf(out, "%i ", getDistance(wantScaffold[j][i - 1], wantScaffold[j][i]));
		fprintf(out, "\n");
	}
	fprintf(out, "Found distance: %i\n", ScaffoldComparer::Compare(foundScaffold, truth));
	fprintf(out, "Want distance: %i\n", ScaffoldComparer::Compare(wantScaffold, truth));

	GraphViz graphWant(wantScaffold, store, want), graphFound(foundScaffold, store, found);
	FILE *wantFile = fopen((config.OutputFileName + "_want.gv").c_str(), "w"), *foundFile = fopen((config.OutputFileName + "_found.gv").c_str(), "w");
	fprintf(wantFile, "%s\n", graphWant.GetString().c_str());
	fprintf(foundFile, "%s\n", graphFound.GetString().c_str());
	fclose(wantFile);
	fclose(foundFile);
}

void testBnB()
{
	vector<bool> tBest(store.ContigCount);
	for (int i = 0; i < store.ContigCount; i++)
		tBest[i] = getOrientation(store[i].Sequence.Name());
	BranchAndBound *bnb = new BranchAndBound(vector<bool>(store.ContigCount, true), tBest, store.ContigCount);
	IterativeSolver *it = new IterativeSolver(vector<bool>(store.ContigCount, true), tBest, store.ContigCount);
	it->Options = bnb->Options = config.Options;
	if (!bnb->Formulate(store) || !bnb->Solve())
	{
		printf("Unable to formulate or solve bnb\n");
		delete bnb;
		delete it;
		return;
	}
	if (!it->Formulate(store) || !it->Solve())
	{
		printf("Unable to formulate or solve itratative\n");
		delete bnb;
		delete it;
		return;
	}
	vector<Scaffold> bnbScaffold = ScaffoldExtractor::Extract(store, *bnb);
	vector<Scaffold> itScaffold = ScaffoldExtractor::Extract(*it);
	vector<Scaffold> truth = ScaffoldExtractor::Extract(store);
	for (int i = 0; i < (int)truth.size(); i++)
		truth[i].NormalizeCoordindates();
	fprintf(out, "Greedy: %.6lf\n", it->GetObjective());
	fprintf(out, "BnB: %.6lf\n", bnb->GetObjective());
	fprintf(out, "Simulation scaffolds: %i\n", (int)truth.size());
	for (int i = 0; i < truth[0].ContigCount(); i++)
		fprintf(out, "%c%i ", (truth[0][i].T ? '-' : '+'), truth[0][i].Id);
	fprintf(out, "\n");
	for (int i = 1; i < truth[0].ContigCount(); i++)
		fprintf(out, "%i ", getDistance(truth[0][i - 1], truth[0][i]));
	fprintf(out, "\n");

	fprintf(out, "Greedy scaffolds: %i, mismatch %i\n", (int)itScaffold.size(), ScaffoldComparer::Compare(itScaffold, truth));
	for (int i = 0; i < itScaffold[0].ContigCount(); i++)
		fprintf(out, "%c%i ", (itScaffold[0][i].T ? '-' : '+'), itScaffold[0][i].Id);
	fprintf(out, "\n");
	for (int i = 1; i < itScaffold[0].ContigCount(); i++)
		fprintf(out, "%i ", getDistance(itScaffold[0][i - 1], itScaffold[0][i]));
	fprintf(out, "\n");

	fprintf(out, "BnB scaffolds: %i, mismatch %i\n", (int)bnbScaffold.size(), ScaffoldComparer::Compare(bnbScaffold, truth));
	for (int i = 0; i < bnbScaffold[0].ContigCount(); i++)
		fprintf(out, "%c%i ", (bnbScaffold[0][i].T ? '-' : '+'), bnbScaffold[0][i].Id);
	fprintf(out, "\n");
	for (int i = 1; i < bnbScaffold[0].ContigCount(); i++)
		fprintf(out, "%i ", getDistance(bnbScaffold[0][i - 1], bnbScaffold[0][i]));
	fprintf(out, "\n");

	fprintf(out, "Distance between scaffolds: %i %i\n", ScaffoldComparer::Compare(itScaffold, bnbScaffold), ScaffoldComparer::Compare(bnbScaffold, itScaffold));

	int count = 0;
	int n = bnb->Incumbent.size();
	for (int i = 0; i < n; i++)
	{
		IloNumVar var = (i < n / 2 ? bnb->alpha[i] : bnb->beta[i - n / 2]);
		double val = bnb->cplex.getValue(var);
		bool sw = ((int)(val + 0.5)) > 0;
		if (bnb->Incumbent[i] != sw)
			count++;
	}
	fprintf(out, "Mismatch: %i Total: %i\n", count, n);
	for (int i = 0; i < n; i++)
	{
		IloNumVar var = (i < n / 2 ? bnb->alpha[i] : bnb->beta[i - n / 2]);
		double val = bnb->cplex.getValue(var);
		bool sw = ((int)(val + 0.5)) > 0;
		if (bnb->Incumbent[i] != sw)
			fprintf(out, "%s %.5lf : %i -> %i\n", (i < n / 2 ? "alpha" : "beta"), bnb->Slack[i], (int)bnb->Incumbent[i], (int)sw);
	}
	delete bnb;
	delete it;
}

void testGA(const GASolver &solver, const string &fileName)
{
	int forwardMismatch = 0, reverseMismatch = 0;
	vector<bool> tBest(store.ContigCount);
	for (int i = 0; i < store.ContigCount; i++)
	{
		tBest[i] = getOrientation(store[i].Sequence.Name());
		if (tBest[i] != solver.T[i])
			forwardMismatch++;
		if (tBest[i] != !solver.T[i])
			reverseMismatch++;
	}
	int mismatch = min(forwardMismatch, reverseMismatch);
	GAIndividual ind(tBest, solver.matrix);
	printf("Found: %.6lf\n", solver.GetObjective());
	printf("Suggested: %.6lf\n", ind.GetObjective());
	printf("Mismatch: %i (%i %i)\n", mismatch, forwardMismatch, reverseMismatch);
	//getFixed(solver.T, tBest);
	vector<Scaffold> single = ScaffoldExtractor::Extract(store);
	for (int i = 0; i < (int)single.size(); i++)
		single[i].NormalizeCoordindates();
	getIterative(solver.T, tBest, single);
	printf("Suggested scaffolds: %i\n", (int)single.size());
	for (int i = 0; i < single[0].ContigCount(); i++)
		printf("%c%i ", (single[0][i].T ? '-' : '+'), single[0][i].Id);
	printf("\n");
	for (int i = 1; i < single[0].ContigCount(); i++)
		printf("%i ", getDistance(single[0][i - 1], single[0][i]));
	printf("\n");
}

void testFormulation()
{
	vector<bool> uBest(store.ContigCount, 1), tBest(store.ContigCount);
	for (int i = 0; i < store.ContigCount; i++)
		tBest[i] = getOrientation(store[i].Sequence.Name());
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

#endif

bool readStore(const string &fileName, DataStore &store)
{
	DataStoreReader reader;
	bool result = reader.Open(fileName) && reader.Read(store);
	reader.Close();
	return result;
}

bool readCoverage(const string &fileName, ReadCoverage &coverage)
{
    ReadCoverageReader reader;
    bool result = reader.Open(fileName) && reader.Read(coverage);
    reader.Close();
    return result;
}

bool outputScaffolds(const string &fileName, const vector<Scaffold> &scaffolds)
{
	FILE *out = fopen(fileName.c_str(), "w");
	if (out == NULL)
		return false;
	int count = scaffolds.size();
	fprintf(out, "%i\n", count);
	for (int j = 0; j < count; j++)
	{
		int size = scaffolds[j].ContigCount();
		for (int i = 0; i < size; i++)
			fprintf(out, "%c%i ", (scaffolds[j][i].T ? '-' : '+'), scaffolds[j][i].Id);
		fprintf(out, "\n");
		for (int i = 0; i < size; i++)
			fprintf(out, "%i ", (int)scaffolds[j][i].X);
		fprintf(out, "\n");
	}
	fclose(out);
	return true;
}

bool outputFastaScaffolds(const string &fileName, const vector<Scaffold> &scaffolds)
{
	FastAWriter writer;
	bool result = writer.Open(fileName) && writer.Write(ScaffoldConverter::ToFasta(store, scaffolds));
	writer.Close();
	return result;
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
		if (config.RemoveAmbiguous)
			cerr << "[+] Removed " << store.RemoveAmbiguous() << " ambiguous links." << endl;
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
		if (config.Erosion > 0)
			cerr << "[i] Erosion removed " << store.Erode(config.Erosion) << " contig links." << endl;
                if (!config.ReadCoverageFileName.empty() && !readCoverage(config.ReadCoverageFileName, coverage))
                    cerr << "[-] Unable to read contig coverage (" << config.ReadCoverageFileName << ")." << endl;
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
		cerr << "[+] Solved the optimization problem." << endl;
		if (!outputFastaScaffolds(config.OutputFileName, ScaffoldExtractor::Extract(solver)))
		{
			cerr << "[-] Unable to output scaffolds (FastA)." << endl;
			return -4;
		}
		if (!config.SolutionOutputFileName.empty() && !outputScaffolds(config.SolutionOutputFileName, ScaffoldExtractor::Extract(solver)))
		{
			cerr << "[-] Unable to output solution." << endl;
			return -5;
		}
		return 0;
	}
	cerr << config.LastError;
	return -1;
}
