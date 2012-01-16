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

#ifndef _GAINDIVIDUAL_H
#define _GAINDIVIDUAL_H
#include "GASolver.h"
#include "GAMatrix.h"
#include <vector>

using namespace std;

class GAIndividual
{
public:
	GAIndividual(int n = 0);
	GAIndividual(const vector<bool> &t, const GAMatrix &matrix);

public:
	double GetObjective() const;
	int GetLength() const;

public:
	void Flip(int i, const GAMatrix &matrix);

public:
	bool operator< (const GAIndividual &other) const;
	bool operator> (const GAIndividual &other) const;
	bool operator== (const GAIndividual &other) const;
	bool operator!= (const GAIndividual &other) const;

private:
	void init(int n);
	void obtainObjectiveValue(const GAMatrix &matrix);
	void initializeGains(const GAMatrix &matrix);
	void updateGains(int k, const GAMatrix &matrix);
	void flip(int k);

public:
	vector<bool> X;
	vector<double> Gain;

private:
	double objectiveValue;
	int length;
};
#endif
