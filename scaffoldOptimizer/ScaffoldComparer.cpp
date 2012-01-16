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

#include "ScaffoldComparer.h"
#include "Helpers.h"
#include <limits>

#include <iostream>

int ScaffoldComparer::Compare(const Scaffold &a, const Scaffold &b)
{
	Scaffold c(b);
	c.Reverse();
	return min(compareOriented(a, b), compareOriented(a, c));
}

int ScaffoldComparer::Compare(const vector<Scaffold> &a, const vector<Scaffold> &b)
{
	int aSize = a.size(), bSize = b.size();
	int mismatch = 0;
	for (int i = 0; i < aSize; i++)
	{
		int minScore = numeric_limits<int>::max();
		for (int j = 0; j < bSize; j++)
			minScore = min(minScore, Compare(a[i], b[j]));
		mismatch += minScore;
	}
	return mismatch;
}

int ScaffoldComparer::compareOriented(const Scaffold &a, const Scaffold &b)
{
	map<int, int> pos;
	int aSize = a.ContigCount(), bSize = b.ContigCount(), mismatch = 0;
	for (int i = 0; i < bSize; i++)
		pos[b[i].Id] = i;
	for (int i = 1; i < aSize; i++)
	{
		ScaffoldContig p = a[i - 1], q = a[i];
		if (pos.find(p.Id) == pos.end() || pos.find(q.Id) == pos.end())
		{
			mismatch++;
			continue;
		}
		if (pos[p.Id] > pos[q.Id])
		{
			mismatch++;
			continue;
		}
		if ((b[pos[p.Id]].T ^ b[pos[q.Id]].T) != (p.T ^ q.T))
		{
			mismatch++;
			continue;
		}
	}
	return mismatch;
}

int ScaffoldComparer::OrientationDistance(const Scaffold &a, const Scaffold &b)
{
	map<int, bool> orientation;
	int aSize = a.ContigCount(), bSize = b.ContigCount(), mismatchForward = 0, mismatchReverse = 0;
	for (int i = 0; i < bSize; i++)
		orientation[b[i].Id] = b[i].T;
	for (int i = 1; i < aSize; i++)
	{
		int id = a[i].Id;
		if (orientation.find(id) == orientation.end())
			mismatchForward++, mismatchReverse++;
		else if (a[i].T == orientation[id])
			mismatchReverse++;
		else
			mismatchForward++;
	}
	return min(mismatchForward, mismatchReverse);
}

int ScaffoldComparer::OrientationDistance(const vector<Scaffold> &a, const vector<Scaffold> &b)
{
	int aSize = a.size(), bSize = b.size();
	int mismatch = 0;
	for (int i = 0; i < aSize; i++)
	{
		int minScore = numeric_limits<int>::max();
		for (int j = 0; j < bSize; j++)
			minScore = min(minScore, OrientationDistance(a[i], b[j]));
		mismatch += minScore;
	}
	return mismatch;
}
