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
	bool PrintMatrix;
	SolverConfiguration Options;
	string InputFileName;
	string OutputFileName;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
