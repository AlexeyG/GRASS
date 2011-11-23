#ifndef _DPGRAPH_H
#define _DPGRAPH_H

#include "DataStore.h"
#include <vector>
#include <set>

using namespace std;

class DPGraph
{
public:
	DPGraph();
	DPGraph(const DataStore &store);
	DPGraph(const DataStore &store, const vector<bool> &t);
	DPGraph(const DataStore &store, const vector<bool> &t, const vector<bool> &l);
	virtual ~DPGraph();

public:
	int FindConnectedComponents(vector< vector<int> > &c);
	int FindBridges(vector< pair<int,int> > &b);
	int RemoveBridge(int a, int b, vector< vector<int> > &c);

public:
	int NVertices;

private:
	void bidirectedConnectedComponentsDFS(int v, vector<int> &colors, int color, vector<int> &component);
	void bidirectedBridgeDFS(int v, int p, vector<bool> &used, vector<int> &tIn, vector<int> &tUp, int &time, vector< pair<int,int> > &b);

private:
	vector< vector<int> > matrix;
	vector< set<int> > list;
};
#endif
