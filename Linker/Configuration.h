#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H
#include <string>
#include <vector>

using namespace std;

class PairedInput
{
public:
	PairedInput(const string &leftFileName, const string &rightFileName, double mean, double std, bool isIllumina, double weight = 1, int mapQ = 0, int minReadLength = 0) : LeftFileName(leftFileName), RightFileName(rightFileName), Mean(mean), Std(std), IsIllumina(isIllumina), Weight(weight), MapQ(mapQ), MinReadLength(minReadLength) {};

public:
	string LeftFileName;
	string RightFileName;
	double Mean;
	double Std;
	bool IsIllumina;
	double Weight;
	int MapQ;
	int MinReadLength;
};

class BWAConfiguration
{
public:
	BWAConfiguration();

public:
	int NumberOfThreads;
	int MaximumHits;
	string IndexCommand;
	string SuffixArrayCommand;
	string AlignSingleEndCommand;
	string TmpPath;
};

class NovoAlignConfiguration
{
public:
	NovoAlignConfiguration();

public:
	string IndexCommand;
	string AlignSingleEndCommand;
	string TmpPath;
};

class SAMToolsConfiguration
{
public:
	SAMToolsConfiguration();

public:
	string ConvertCommand;
	string TmpPath;
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
	int MaximumLinkHits;
	BWAConfiguration BWAConfig;
	NovoAlignConfiguration NovoAlignConfig;
	SAMToolsConfiguration SAMToolsConfig;
	vector<PairedInput> PairedReadInputs;
	string LastError;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
