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

#include "Globals.h"
#include "GraphViz.h"
#include "Helpers.h"
#include "IterativeSolver.h"
#include "ExtendedFixedMIQPSolver.h"
#include <limits>

GraphViz::GraphViz(const vector<Scaffold> &scaffold, const DataStore &store, const IterativeSolver &solver)
	: scaffold(scaffold), T(store.ContigCount, false), pos(store.ContigCount)
{
	header();
	for (int i = 0; i < (int)scaffold.size(); i++)
	{
		outputScaffold(scaffold[i], i);
		for (int j = 0; j < scaffold[i].ContigCount(); j++)
		{
			T[scaffold[i][j].Id] = scaffold[i][j].T;
			pos[scaffold[i][j].Id] = j;
		}
	}
	outputLinks(store, solver);
	footer();
}

string GraphViz::GetString()
{
	return out.str();
}

void GraphViz::putline(const char *format, ...)
{
	char buf[MaxLine];
	va_list list;
	va_start(list, format);
	vsprintf(buf, format, list);
	va_end(list);
	out << string(buf) << endl;
}

void GraphViz::header()
{
	putline("digraph scaffold");
	putline("{");
	putline("\trankdir=LR;");
	putline("\tnode [shape=house, height=0.3, fixedsize=true];");
}

void GraphViz::footer()
{
	putline("}");
}

void GraphViz::outputScaffold(const Scaffold &scaffold, int id)
{
	int size = scaffold.ContigCount();
	putline("\tsubgraph scaffold_%i", id);
	putline("\t{");
	putline("\t\tlabel = \"Scaffold %i\";", id + 1);
	for (int i = 0; i < size; i++)
		putline("\t\t%i [orientation=%i];", scaffold[i].Id + 1, (scaffold[i].T ? 270 : 90));
	putline("\t\tedge [style=invis];");
	for (int i = 1; i < size; i++)
		putline("\t\t%i->%i;", scaffold[i - 1].Id + 1, scaffold[i].Id + 1);
	putline("\t}");
}

void GraphViz::outputLinks(const DataStore &store, const IterativeSolver &solver)
{
	double maxWeight = -Helpers::Inf;
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
	{
		int a = it->first.first, b = it->first.second;
		if ((T[a] ^ T[b]) != it->second.EqualOrientation)
			if (maxWeight < it->second.Weight)
				maxWeight = it->second.Weight;
	}
	int num = 0;
	putline("\tedge [style=solid, constraint=false, fontsize=10];");
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
	{
		int a = it->first.first, b = it->first.second;
		if ((T[a] ^ T[b]) != it->second.EqualOrientation)
		{
			double xi = solver.GetDistanceSlack(num), delta = solver.GetOrderSlack(num);
			bool dashed = false;
			if (xi > ExtendedFixedMIQPSolver::DesiredDistanceSlackMax - Helpers::Eps || delta > ExtendedFixedMIQPSolver::DesiredOrderSlackMax - Helpers::Eps)
				dashed = true;

			string colorXi = getColor(xi, ExtendedFixedMIQPSolver::DesiredDistanceSlackMax), colorDelta = getColor(delta, ExtendedFixedMIQPSolver::DesiredOrderSlackMax);

			if ((!T[a] && it->second.ForwardOrder) || (T[a] && !it->second.ForwardOrder))
				putline("\t%i->%i [label=\"%.2lf&plusmn;%.2lf\", penwidth=%.5lf, style=%s, color=\"%s:white:%s\"];", a + 1, b + 1, it->second.Mean, it->second.Std, it->second.Weight / maxWeight * (MaxPenWidth - MinPenWidth) + MinPenWidth, (dashed ? "dashed" : "solid"), colorDelta.c_str(), colorXi.c_str());
			else
				putline("\t%i->%i [label=\"%.2lf&plusmn;%.2lf\", penwidth=%.5lf, style=%s, color=\"%s:white:%s\"];", b + 1, a + 1, it->second.Mean, it->second.Std, it->second.Weight / maxWeight * (MaxPenWidth - MinPenWidth) + MinPenWidth, (dashed ? "dashed" : "solid"), colorDelta.c_str(), colorXi.c_str());
			num++;
		}
	}
}

string GraphViz::getColor(double s, double max)
{
	char buf[MaxLine];
	if (s < Helpers::Eps)
		return "black";
	if (s > max - Helpers::Eps)
		return "red";
	unsigned char redness = (s * 255) / max;
	sprintf(buf, "#%2X0000", redness);
	return string(buf);
}
