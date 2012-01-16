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

#ifndef _GAMATRIX_H
#define _GAMATRIX_H
#include <vector>

using namespace std;

class GAMatrix
{
public:
	typedef vector< vector<double> > Matrix;

public:
	GAMatrix(int n = 0);

public:
	void CalculatePositions();

public:
	vector<double> &operator[] (int i);
	const vector<double> &operator[] (int i) const;

private:
	int size;
	Matrix matrix;

public:
	vector< vector<int> > Pos;
};
#endif
