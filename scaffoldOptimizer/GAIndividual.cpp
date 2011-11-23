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
