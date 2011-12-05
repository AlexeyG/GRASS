#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H
#include <string>
#include <vector>
#include "AlignerConfiguration.h"

using namespace std;

class Configuration
{
public:
	Configuration();

public:
	bool ProcessCommandLine(int argc, char *argv[]);

public:
	bool Success;
	string ScaffoldFileName;
	string ReferenceFileName;
        int MinBases;
        MummerConfiguration MummerConfig;
	string LastError;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
