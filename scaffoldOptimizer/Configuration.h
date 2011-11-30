#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H
#include <string>
#include <vector>
#include <sstream>
#include "SolverConfiguration.h"

using namespace std;

class Configuration
{
public:
	Configuration();

public:
	bool ProcessCommandLine(int argc, char *argv[]);

public:
	bool Success;
	string LastError;
	bool RemoveAmbiguous;
	bool Sort;
	bool Bundle;
	bool BundlePerGroup;
	bool BundleAmbiguous;
	double BundleDistance;
	double Erosion;
        double ExpectedCoverage;
        double UniquenessFCutoff;
	bool PrintMatrix;
	SolverConfiguration Options;
	string InputFileName;
        string ReadCoverageFileName;
	string OutputFileName;
	string SolutionOutputFileName;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
