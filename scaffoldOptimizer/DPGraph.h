/*
 * scaffoldOptimizer : solves the MIQP optimization and produces linear scaffold
 * sequences.
 * Copyright (C) 2011  Alexey Gritsenko
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 * 
 * 
 * Email: a.gritsenko@tudelft.nl
 * Mail: Delft University of Technology
 *       Faculty of Electrical Engineering, Mathematics, and Computer Science
 *       Department of Mediamatics
 *       P.O. Box 5031
 *       2600 GA, Delft, The Netherlands
 */


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
