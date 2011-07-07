#include "GAMatrix.h"

GAMatrix::GAMatrix(int n)
	: size(n), matrix(n, vector<double>(n, 0)), Pos(n, vector<int>())
{
}

void GAMatrix::CalculatePositions()
{
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) // check that we don't need just a half of it
			if (matrix[i][j] != 0)
				Pos[i].push_back(j);
}

vector<double> &GAMatrix::operator[] (int i)
{
	return matrix[i];
}

const vector<double> &GAMatrix::operator[] (int i) const
{
	return matrix[i];
}
