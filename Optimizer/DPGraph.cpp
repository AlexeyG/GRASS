#include "DPGraph.h"

DPGraph::DPGraph()
{
}

DPGraph::DPGraph(const DataStore &store)
{
	NVertices = store.ContigCount;
	for (int i = 0; i < NVertices; i++)
		for (int j = 0; j < NVertices; j++)
		{
			DataStore::LinkRange range = store(i, j);
			for (DataStore::LinkMap::const_iterator it = range.first; it != range.second; it++)
				this->Update(i, j, (*this)(i, j) + 1);
		}
}

DPGraph::~DPGraph()
{
}

int DPGraph::operator() (int i, int j)
{
	pair<int, int> q(i, j);
	Matrix::iterator found = data.find(q);
	if (found != data.end())
		return found->second;
	else
		return 0;
}

int DPGraph::Update(int i, int j, int newValue)
{
	pair<int, int> q(i, j);
	if (newValue == 0)
		data.erase(q);
	else
		data[q] = newValue;
	return newValue;
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
	int ab = (*this)(a, b), ba = (*this)(b, a);
	Update(a, b, 0); Update(b, a, 0);
	int res = FindConnectedComponents(c);
	Update(a, b, ab); Update(b, a, ba);
	return res;
}

void DPGraph::bidirectedConnectedComponentsDFS(int v, vector<int> &colors, int color, vector<int> &component)
{
	component.push_back(v);
	colors[v] = color;
	for (int i = 0; i < NVertices; i++)
		if (((*this)(v, i) > 0 || (*this)(i, v) > 0) && colors[i] < 0)
			bidirectedConnectedComponentsDFS(i, colors, color, component);
}

void DPGraph::bidirectedBridgeDFS(int v, int p, vector<bool> &used, vector<int> &tIn, vector<int> &tUp, int &time, vector< pair<int,int> > &b)
{
	used[v] = true;
	tIn[v] = tUp[v] = time++;
	for (int i = 0; i < NVertices; i++)
	{
		if (i == p) continue;
		if ((*this)(v, i) > 0 || (*this)(i, v) > 0)
		{
			if (used[i])
				tUp[v] = min(tUp[v], tIn[i]);
			else
			{
				bidirectedBridgeDFS(i, v, used, tIn, tUp, time, b);
				tUp[v] = min(tUp[v], tUp[i]);
				if (tUp[i] > tIn[v])
					b.push_back(pair<int,int>(v, i));
			}
		}
	}
}
