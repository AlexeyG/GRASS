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
