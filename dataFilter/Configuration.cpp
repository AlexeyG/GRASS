#include "Configuration.h"
#include "Defines.h"
#include "Helpers.h"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <sstream>

using namespace std;

PairedInput::PairedInput(const string &leftFileName, const string &rightFileName, const string &outputPrefix)
	: LeftFileName(leftFileName), RightFileName(rightFileName), OutputPrefix(outputPrefix)
{
}

// Construtor with default configuration parameter settings.
Configuration::Configuration()
{
	Success = false;
	InputFileName = "";
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
			else if (!strcmp("-paired", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -paired: must have three arguments." << endl;
					this->Success = false;
					break;
				}
				i++; string leftFileName = argv[i];
				i++; string rightFileName = argv[i];
				i++; string outputPrefix = argv[i];
				PairedFilter.push_back(PairedInput(leftFileName, rightFileName, outputPrefix));
			}
			else if (i == argc - 1)
				this->InputFileName = argv[argc - 1];
			else
			{
				serr << "[-] Unknown argument: " << argv[i] << endl;
				this->Success = false;
				break;
			}
			i++;
		}
		if (this->InputFileName == "")
		{
			serr << "[-] No input file specified." << endl;
			this->Success = false;
		}
	}
	if (!this->Success)
		LastError = serr.str();
	return this->Success;
}

void Configuration::printHelpMessage(stringstream &serr)
{
	serr << "[i] Scaffold optimizer version " << VERSION << " (" << DATE << ")" << endl;
	serr << "[i] By " << AUTHOR << endl;
	serr << "[i] Usage: scaffoldOptimizer [arguments] <scaffold.opt>" << endl;
	serr << "[i] -help                                               Print this message and exit." << endl;
	serr << "[i] -paired <left-file> <right-file> <output-prefix>    Filter paired reads and output them with given prefix." << endl;
	serr << "[i] -output [output filename]                           Output filename for optimzation information. [scaffold.fasta]" << endl;
}
