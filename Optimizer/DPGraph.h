#ifndef _DPGRAPH_H
#define _DPGRAPH_H

#include "DataStore.h"
#include <map>
#include <vector>

using namespace std;

class DPGraph
{
public:
	DPGraph();
	DPGraph(const DataStore &store);
	virtual ~DPGraph();

public:
	int operator() (int i, int j);
	int Update(int i, int j, int newValue);
	int FindConnectedComponents(vector< vector<int> > &c);
	int FindBridges(vector< pair<int,int> > &b);
	int RemoveBridge(int a, int b, vector< vector<int> > &c);

public:
	int NVertices;

private:
	void bidirectedConnectedComponentsDFS(int v, vector<int> &colors, int color, vector<int> &component);
	void bidirectedBridgeDFS(int v, int p, vector<bool> &used, vector<int> &tIn, vector<int> &tUp, int &time, vector< pair<int,int> > &b);

private:
	typedef map<pair<int, int>, int> Matrix;

private:
	map<pair<int, int>, int> data;
};
#endif
