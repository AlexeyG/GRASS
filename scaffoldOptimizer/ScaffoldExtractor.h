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
