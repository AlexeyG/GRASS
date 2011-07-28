#ifndef _SCAFFOLDEXTRACTOR_H
#define _SCAFFOLDEXTRACTOR_H

#include <vector>
#include <algorithm>

using namespace std;

class DataStore;
class IterativeSolver;
class GASolver;
class FixedMIQPSolver;
class BranchAndBound;
class EMSolver;
class DPSolver;
class Contig;

struct ScaffoldContig
{
public:
	ScaffoldContig(int id, bool t, double x, int len);

public:
	int Id;
	bool T;
	double X;
	int Length;

public:
	bool operator< (const ScaffoldContig &b) const;
};

class Scaffold
{
public:
	Scaffold();

public:
	void AddContig(int id, bool t, double x, int len);
	void AddContig(const ScaffoldContig &contig);
	void Sort();
	void NormalizeCoordindates();
	void Reverse();
	void ApplyTransform(const vector<int> &transform);
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
	static vector<Scaffold> Extract(const DataStore &store, const FixedMIQPSolver &solver);
	static vector<Scaffold> Extract(const DataStore &store, const BranchAndBound &solver);
	static vector<Scaffold> Extract(const EMSolver &solver);
	static vector<Scaffold> Extract(const DPSolver &solver);

private:
	static bool getOrientation(const Contig &contig, bool &orientation, double &position);
	static void extractSingleScaffold(const DataStore &store, vector<Scaffold> &ans);
	static void extractOrientedScaffold(const DataStore &store, vector<Scaffold> &ans);
};
#endif
