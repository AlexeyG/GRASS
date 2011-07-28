#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H
#include <string>
#include <vector>
#include <sstream>

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
	int Kmer;
	bool Verbose;
	string InputFileName;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
