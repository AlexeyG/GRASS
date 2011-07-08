#ifndef _SCAFFOLDEXTRACTOR_H
#define _SCAFFOLDEXTRACTOR_H

#include <vector>
#include <algorithm>
#include "DataStore.h"
#include "IterativeSolver.h"
#include "FixedMIQPSolver.h"
#include "GASolver.h"

using namespace std;

struct ScaffoldContig
{
public:
	ScaffoldContig(int id, bool t, double x);

public:
	int Id;
	bool T;
	double X;

public:
	bool operator< (const ScaffoldContig &b) const;
};

class Scaffold
{
public:
	Scaffold();

public:
	void AddContig(int id, bool t, double x);
	void AddContig(const ScaffoldContig &contig);
	void Sort();
	const ScaffoldContig &operator[] (int i) const;
	int ContigCount() const;

private:
	vector<ScaffoldContig> contigs;
};

class ScaffoldExtractor
{
public:
	static vector<Scaffold> Extract(const DataStore &store, bool single = true);
	static vector<Scaffold> Extract(const IterativeSolver &sovler);
	static vector<Scaffold> Extract(const DataStore &store, const GASolver &solver);
	static vector<Scaffold> Extract(const DataStore &store, const FixedMIQPSolver &sovler);

private:
	static bool getOrientation(const Contig &contig, bool &orientation, double &position);
	static void extractSingleScaffold(const DataStore &store, vector<Scaffold> &ans);
	static void extractOrientedScaffold(const DataStore &store, vector<Scaffold> &ans);
};
#endif
