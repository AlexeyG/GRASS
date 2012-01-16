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

#ifndef _RANDOMIZEDGREEDYINITIALIZER_H
#define _RANDOMIZEDGREEDYINITIALIZER_H
#include "GAIndividual.h"
#include "GAMatrix.h"
#include <vector>

using namespace std;

class RandomizedGreedyInitializer
{
public:
	RandomizedGreedyInitializer(int n, const GAMatrix &matrix);
	
public:
	GAIndividual MakeSolution(const GAMatrix &matrix);

private:
	void initializeGains(const GAMatrix &matrix);
	void updateGains(int k, bool value, const GAMatrix &matrix);
	void updateList(int k, bool value);
	void flip(int k, bool value);

	/*bool checkGains(const GAMatrix &matrix) const;
	void print() const;*/

private:
	int length, unset;
	vector<double> x;
	vector<double> gainZero, gainOne;
	vector<bool> selected;
};
#endif
