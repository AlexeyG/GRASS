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

#include "GAIndividual.h"
#include "Helpers.h"

GAIndividual::GAIndividual(int n)
{
	init(n);
}

GAIndividual::GAIndividual(const vector<bool> &t, const GAMatrix &matrix)
{
	init(t.size());
	for (int i = 0; i < length; i++)
		X[i] = t[i];
	obtainObjectiveValue(matrix);
	initializeGains(matrix);
}

double GAIndividual::GetObjective() const
{
	return objectiveValue;
}

int GAIndividual::GetLength() const
{
	return length;
}

void GAIndividual::Flip(int i, const GAMatrix &matrix)
{
	objectiveValue += Gain[i];
	updateGains(i, matrix);
	flip(i);
}

bool GAIndividual::operator< (const GAIndividual &other) const
{
	return objectiveValue < other.objectiveValue;
}

bool GAIndividual::operator> (const GAIndividual &other) const
{
	return objectiveValue > other.objectiveValue;
}

bool GAIndividual::operator== (const GAIndividual &other) const
{
	return fabs(objectiveValue - other.objectiveValue) < Helpers::Eps;
}

bool GAIndividual::operator!= (const GAIndividual &other) const
{
	return fabs(objectiveValue - other.objectiveValue) >= Helpers::Eps;
}

void GAIndividual::init(int n)
{
	objectiveValue = 0;
	length = n;
	X.resize(n);
	Gain.resize(n);
}

void GAIndividual::obtainObjectiveValue(const GAMatrix &matrix)
{
	objectiveValue = matrix[length][length]; // constant
	for (int i = 0; i < length; i++)
	{
		objectiveValue += matrix[i][i] * (int)X[i];
		for (vector<int>::const_iterator j = matrix.Pos[i].begin(); j != matrix.Pos[i].end(); j++)
			if (*j > i)
				objectiveValue += 2 * matrix[i][*j] * (int)X[i] * (int)X[*j];
	}
}

void GAIndividual::initializeGains(const GAMatrix &matrix)
{
	for (int i = 0; i < length; i++)
	{
		double sum = matrix[i][i] * ((int)(!X[i]) - (int)X[i]);
		for (vector<int>::const_iterator j = matrix.Pos[i].begin(); j != matrix.Pos[i].end(); j++)
			if (*j != i)
				sum += 2 * matrix[*j][i] * (int)X[*j] * ((int)(!X[i]) - (int)X[i]);
		Gain[i] = sum;
	}
}

void GAIndividual::updateGains(int k, const GAMatrix &matrix)
{
	Gain[k] = -Gain[k];
	for (vector<int>::const_iterator i = matrix.Pos[k].begin(); i != matrix.Pos[k].end(); i++)
		if (*i != k)
			Gain[*i] += 2 * matrix[*i][k] * ((int)!X[*i] - (int)X[*i]) * ((int)!X[k] - (int)X[k]);
}

void GAIndividual::flip(int k)
{
	X[k].flip();
}
