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

#include "RandomizedGreedyInitializer.h"
#include "Helpers.h"
#include <cstdlib>

RandomizedGreedyInitializer::RandomizedGreedyInitializer(int n, const GAMatrix &matrix)
	: length(n), unset(n), x(n, 0.5), gainZero(n), gainOne(n), selected(n, false)
{
	initializeGains(matrix);
}

GAIndividual RandomizedGreedyInitializer::MakeSolution(const GAMatrix &matrix)
{
	if (unset > 0)
	{
		int k = rand() % length;
		bool l = rand() < RAND_MAX / 2;
		updateGains(k, l, matrix);
		updateList(k, l);
		flip(k, l);
		while (unset > 0)
		{
			int k0 = -1;
			int k1 = -1;
			for (int i = 0; i < length; i++)
				if (!selected[i])
				{
					if (k0 < 0 || gainZero[i] > gainZero[k0])
						k0 = i;
					if (k1 < 0 || gainOne[i] > gainOne[k1])
						k1 = i;
				}
			double sum = gainZero[k0] + gainOne[k1];
			double p = (sum < Helpers::Eps ? 0.5 : gainZero[k0] / sum);
			if (rand() < (int)(RAND_MAX * p))
				k = k0, l = false;
			else
				k = k1, l = true;

			updateGains(k, l, matrix);
			updateList(k, l);
			flip(k, l);
		}
	}
	vector<bool> t(length);
	for (int i = 0; i < length; i++)
		t[i] = x[i] == 1;
	return GAIndividual(t, matrix);
}

void RandomizedGreedyInitializer::initializeGains(const GAMatrix &matrix)
{
	for (int i = 0; i < length; i++)
	{
		gainZero[i] = -0.25 * matrix[i][i];
		gainOne[i] = 0.75 * matrix[i][i];
		for (vector<int>::const_iterator j = matrix.Pos[i].begin(); j != matrix.Pos[i].end(); j++)
			if (*j != i)
			{
				gainZero[i] -= matrix[i][*j] * x[*j];
				gainOne[i] += matrix[i][*j] * x[*j];
			}
	}
}

void RandomizedGreedyInitializer::updateGains(int k, bool value, const GAMatrix &matrix)
{
	if (value)
	{
		for (vector<int>::const_iterator i = matrix.Pos[k].begin(); i != matrix.Pos[k].end(); i++)
			if (*i != k)
			{
				gainZero[*i] -= 0.5 * matrix[*i][k];
				gainOne[*i] += 0.5 * matrix[*i][k];
			}
	}
	else
		for (vector<int>::const_iterator i = matrix.Pos[k].begin(); i != matrix.Pos[k].end(); i++)
			if (*i != k)
			{
				gainZero[*i] += 0.5 * matrix[*i][k];
				gainOne[*i] -= 0.5 * matrix[*i][k];
			}
}

void RandomizedGreedyInitializer::updateList(int k, bool value)
{
	unset--;
	selected[k] = true;
}

void RandomizedGreedyInitializer::flip(int k, bool value)
{
	if (value)
		x[k] = 1;
	else
		x[k] = 0;
}
