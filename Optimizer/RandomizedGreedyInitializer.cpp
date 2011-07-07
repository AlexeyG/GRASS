#include "RandomizedGreedyInitializer.h"
#include "Helpers.h"

RandomizedGreedyInitializer::RandomizedGreedyInitializer(int n, const GAMatrix &matrix)
	: length(n), unset(n), x(n, 0.5), gainZero(n), gainOne(n), selected(n, false)
{
	/*for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < length; j++)
			printf("%7.3lf ", matrix[i][j]);
		cout << endl;
	}
	cout << "Got assignment" << endl;*/
	initializeGains(matrix);
	//cout << "Got init gains" << endl;
	//print();
	//if (!checkGains(matrix))
//		throw std::exception();
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
		//cout << "Flip " << k << " to " << l << endl;
		//print();
		//if (!checkGains(matrix))
		//	throw std::exception();
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
			//cout << "Flip " << k << " to " << l << endl;
			//print();
			//if (!checkGains(matrix))
			//	throw std::exception();
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

/*bool RandomizedGreedyInitializer::checkGains(const GAMatrix &matrix) const
{
	for (int j = 0; j < length; j++)
	{
		if (selected[j])
			continue;
		vector<double> y(x);
		y[j] = 0;
		double gain = 0;
		for (int l = 0; l < length; l++)
			for (int m = 0; m < length; m++)
				gain += matrix[l][m] * (y[l] * y[m] - x[l] * x[m]);
		if (fabs(gain - gainZero[j]) > Helpers::Eps)
		{
			cout << "Error at zero-gain: " << j << endl;
			cout << gainZero[j] << " vs. " << gain << endl;
			return false;
		}
		else
			cout << "Zero-gain " << j << " is OK" << endl;
	}

	for (int j = 0; j < length; j++)
	{
		if (selected[j])
			continue;
		vector<double> y(x);
		y[j] = 1;
		double gain = 0;
		for (int l = 0; l < length; l++)
			for (int m = 0; m < length; m++)
				gain += matrix[l][m] * (y[l] * y[m] - x[l] * x[m]);
		if (fabs(gain - gainOne[j]) > Helpers::Eps)
		{
			cout << "Error at one-gain: " << j << endl;
			cout << gainOne[j] << " vs. " << gain << endl;
			return false;
		}
		else
			cout << "One-gain " << j << " is OK" << endl;
	}
	return true;
}

void RandomizedGreedyInitializer::print() const
{
	for (int i = 0; i < length; i++)
		printf("%5.1lf ", x[i]);
	cout << endl;
	cout << "Zero gains:" << endl;
	for (int i = 0; i < length; i++)
		printf("%5.2lf ", gainZero[i]);
	cout << endl;
	cout << "One gains:" << endl;
	for (int i = 0; i < length; i++)
		printf("%5.2lf ", gainOne[i]);
	cout << endl;
}*/

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
