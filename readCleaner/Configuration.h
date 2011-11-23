#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H
#include <string>
#include <vector>
#include "AlignerConfiguration.h"

using namespace std;

class PairedInput
{
public:
	PairedInput(const string &leftFileName, const string &rightFileName, const string &outputPrefix, bool isIllumina) : LeftFileName(leftFileName), RightFileName(rightFileName), OutputPrefix(outputPrefix), IsIllumina(isIllumina) {};

public:
	string LeftFileName;
	string RightFileName;
	string OutputPrefix;
	bool IsIllumina;
};

class Configuration
{
public:
	Configuration();

public:
	bool ProcessCommandLine(int argc, char *argv[]);

public:
	bool Success;
	string ReferenceFileName;
	BWAConfiguration BWAConfig;
	NovoAlignConfiguration NovoAlignConfig;
	SAMToolsConfiguration SAMToolsConfig;
	vector<PairedInput> PairedReadInputs;
	string LastError;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
