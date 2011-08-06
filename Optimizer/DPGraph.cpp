#include "DPGraph.h"

DPGraph::DPGraph()
{
}

DPGraph::DPGraph(const DataStore &store)
{
	NVertices = store.ContigCount;
	matrix.assign(NVertices, vector<int>(NVertices, 0));
	list.assign(NVertices, set<int>());
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
	{
		int i = it->first.first, j = it->first.second;
		matrix[i][j]++, matrix[j][i]++;
		list[i].insert(j), list[j].insert(i);
	}
}

DPGraph::DPGraph(const DataStore &store, const vector<bool> &t)
{
	NVertices = store.ContigCount;
	matrix.assign(NVertices, vector<int>(NVertices, 0));
	list.assign(NVertices, set<int>());
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
	{
		int i = it->first.first, j = it->first.second;
		if ((t[i] ^ t[j]) == it->second.EqualOrientation)
			continue;
		matrix[i][j]++, matrix[j][i]++;
		list[i].insert(j), list[j].insert(i);
	}
}

DPGraph::DPGraph(const DataStore &store, const vector<bool> &t, const vector<bool> &l)
{
	NVertices = store.ContigCount;
	matrix.assign(NVertices, vector<int>(NVertices, 0));
	list.assign(NVertices, set<int>());
	int num = 0;
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
	{
		int i = it->first.first, j = it->first.second;
		if ((t[i] ^ t[j]) == it->second.EqualOrientation)
			continue;
		if (l[num++])
		{
			matrix[i][j]++, matrix[j][i]++;
			list[i].insert(j), list[j].insert(i);
		}
	}
}

DPGraph::~DPGraph()
{
}

int DPGraph::FindConnectedComponents(vector< vector<int> > &c)
{
	vector<int> colors(NVertices, -1);
	int color = 0;
	for (int i = 0; i < NVertices; i++)
		if (colors[i] == -1)
		{
			c.push_back(vector<int>());
			bidirectedConnectedComponentsDFS(i, colors, color, c[color]);
			color++;
		}
	return color;
}

int DPGraph::FindBridges(vector< pair<int,int> > &b)
{
	vector<bool> used(NVertices, false);
	vector<int> tIn(NVertices, 0);
	vector<int> tUp(NVertices, 0);
	int time = 0;
	for (int i = 0; i < NVertices; i++)
		if (!used[i]) bidirectedBridgeDFS(i, -1, used, tIn, tUp, time, b);
	return b.size();
}

int DPGraph::RemoveBridge(int a, int b, vector< vector<int> > &c)
{
	int ab = matrix[a][b];
	matrix[a][b] = matrix[b][a] = 0;
	int res = FindConnectedComponents(c);
	matrix[a][b] = matrix[b][a] = ab;
	return res;
}

void DPGraph::bidirectedConnectedComponentsDFS(int v, vector<int> &colors, int color, vector<int> &component)
{
	component.push_back(v);
	colors[v] = color;
	for (set<int>::const_iterator i = list[v].begin(); i != list[v].end(); i++)
		if (matrix[v][*i] > 0 && colors[*i] < 0)
			bidirectedConnectedComponentsDFS(*i, colors, color, component);
}

void DPGraph::bidirectedBridgeDFS(int v, int p, vector<bool> &used, vector<int> &tIn, vector<int> &tUp, int &time, vector< pair<int,int> > &b)
{
	used[v] = true;
	tIn[v] = tUp[v] = time++;
	for (set<int>::const_iterator i = list[v].begin(); i != list[v].end(); i++)
	{
		if (*i == p) continue;
		if (matrix[v][*i] > 0)
		{
			if (used[*i])
				tUp[v] = min(tUp[v], tIn[*i]);
			else
			{
				bidirectedBridgeDFS(*i, v, used, tIn, tUp, time, b);
				tUp[v] = min(tUp[v], tUp[*i]);
				if (tUp[*i] > tIn[v])
					b.push_back(pair<int,int>(v, *i));
			}
		}
	}
}
