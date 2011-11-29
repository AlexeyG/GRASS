#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H
#include <string>
#include <vector>
#include "AlignerConfiguration.h"

using namespace std;

class PairedInput
{
public:
	PairedInput(const string &leftFileName, const string &rightFileName, double mean, double std, bool isIllumina, double weight = 1, int mapQ = 0, int minReadLength = 0, int maxEditDistance = 10000) : LeftFileName(leftFileName), RightFileName(rightFileName), Mean(mean), Std(std), IsIllumina(isIllumina), Weight(weight), MapQ(mapQ), MinReadLength(minReadLength), MaxEditDistance(maxEditDistance) {};

public:
	string LeftFileName;
	string RightFileName;
	double Mean;
	double Std;
	bool IsIllumina;
	double Weight;
	int MapQ;
	int MinReadLength;
	int MaxEditDistance;
};

class Configuration
{
public:
	Configuration();

public:
	bool ProcessCommandLine(int argc, char *argv[]);

public:
	bool Success;
	string InputFileName;
	string OutputFileName;
	string ReadCoverageFileName;
        int MaximumLinkHits;
	double NoOverlapDeviation;
	BWAConfiguration BWAConfig;
	NovoAlignConfiguration NovoAlignConfig;
	SAMToolsConfiguration SAMToolsConfig;
	vector<PairedInput> PairedReadInputs;
	string LastError;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
