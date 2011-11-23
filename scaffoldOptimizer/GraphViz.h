#ifndef _GRAPHVIZ_H
#define _GRAPHVIZ_H
#include "ScaffoldExtractor.h"
#include "DataStore.h"
#include <vector>
#include <string>
#include <sstream>

using namespace std;

class IterativeSolver;

class GraphViz
{
public:
	GraphViz(const vector<Scaffold> &scaffold, const DataStore &store, const IterativeSolver &solver);

public:
	string GetString();

public:
	static const double MaxPenWidth = 2;
	static const double MinPenWidth = 0.1;

private:
	void putline(const char *format, ...);
	void header();
	void footer();
	void outputScaffold(const Scaffold &scaffold, int id);
	void outputLinks(const DataStore &store, const IterativeSolver &solver);
	static string getColor(double s, double max);

protected:
	vector<Scaffold> scaffold;
	vector<bool> T;
	vector<int> pos;

private:
	ostringstream out;
};
#endif
