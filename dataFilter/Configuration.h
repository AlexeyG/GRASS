#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class PairedInput
{
public:
	PairedInput(const string &leftFileName, const string &rightFileName, const string &outputPrefix);

public:
	string LeftFileName;
	string RightFileName;
	string OutputPrefix;
};

class Configuration
{
public:
	Configuration();

public:
	bool ProcessCommandLine(int argc, char *argv[]);

public:
	bool Success;
	string LastError;
	vector<PairedInput> PairedFilter;
	string InputFileName;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
