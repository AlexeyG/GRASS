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
