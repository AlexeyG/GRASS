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
