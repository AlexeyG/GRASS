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
