/*
 * dataSimulator : a support tool used in debugging to create simulated scaffolding
 * problems given a reference sequence and paired reads.
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

#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H
#include <string>
#include <vector>

using namespace std;

class Split
{
public:
	Split(int mean, int variance, int count = 1) : Mean(mean), Variance(variance), Count(count) {};

public:
	int Mean;
	int Variance;
	int Count;
};

class Configuration
{
public:
	Configuration();

public:
	bool Success;
	bool AllowContigOverlap;
	bool FlipOrientation;
	bool ShuffleContigs;
	int Limit;
	string InputFileName;
	string OutputFileName;
	vector<Split> Splits;
};

#endif
