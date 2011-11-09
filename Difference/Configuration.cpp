#include "Configuration.h"
#include "Defines.h"
#include "Helpers.h"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <sstream>

using namespace std;

// Construtor with default configuration parameter settings.
Configuration::Configuration()
{
	Success = false;
	AFileName = "";
	BFileName = "";
	CFileName = "";
}

// Parses command line arguments. Returns true if successful.
bool Configuration::ProcessCommandLine(int argc, char *argv[])
{
	this->Success = true;
	stringstream serr;

	if (argc == 1)
	{
		serr << "[-] Not enough arguments. Consult -help." << endl;
		this->Success = false;
	}
	else
	{
		int i = 1;
		while (i < argc)
		{
			if (!strcmp("-help", argv[i]) || !strcmp("-h", argv[i]))
			{
				printHelpMessage(serr);
				this->Success = false;
				break;
			}
			else if (i == argc - 3)
				this->AFileName = argv[argc - 3];
			else if (i == argc - 2)
				this->BFileName = argv[argc - 2];
			else if (i == argc - 1)
				this->CFileName = argv[argc - 1];
			else
			{
				serr << "[-] Unknown argument: " << argv[i] << endl;
				this->Success = false;
				break;
			}
			i++;
		}
		if (this->AFileName == "")
		{
			serr << "[-] No first input file specified." << endl;
			this->Success = false;
		}
		if (this->BFileName == "")
		{
			serr << "[-] No second input file specified." << endl;
			this->Success = false;
		}
		if (this->CFileName == "")
		{
			serr << "[-] No output file specified." << endl;
			this->Success = false;
		}
	}
	if (!this->Success)
		LastError = serr.str();
	return this->Success;
}

void Configuration::printHelpMessage(stringstream &serr)
{
	serr << "[i] Read set difference calculator out = (A - B). Version " << VERSION << " (" << DATE << ")" << endl;
	serr << "[i] By " << AUTHOR << endl;
	serr << "[i] Usage: readDiff [arguments] <A.fasta> <B.fasta> <out.fasta>" << endl;
	serr << "[i] -help                                               Print this message and exit." << endl;
}
