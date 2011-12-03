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
	CoverageFileName = "";
        ContigFileName = "";
	DepthFileName = "";
}

// Parses command line arguments. Returns true if successful.
bool Configuration::ProcessCommandLine(int argc, char *argv[])
{
	this->Success = true;
	stringstream serr;

	if (argc < 3)
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
                        else if (!strcmp("-d", argv[i]) || !strcmp("-depth", argv[i]))
                        {
                            if (argc - i - 1 < 1)
                            {
                                serr << "[-] Parsing error in -depth: must have an argument." << endl;
                                this->Success = false;
                                break;
                            }
                            i++;
                            DepthFileName = argv[i];
                        }
                        else if (i == argc - 2)
				this->ContigFileName = argv[argc - 2];
			else if (i == argc - 1)
				this->CoverageFileName = argv[argc - 1];
			else
			{
				serr << "[-] Unknown argument: " << argv[i] << endl;
				this->Success = false;
				break;
			}
			i++;
		}
		if (this->ContigFileName == "")
		{
			serr << "[-] No contig input file specified." << endl;
			this->Success = false;
		}
                if (this->CoverageFileName == "")
		{
			serr << "[-] No coverage input file specified." << endl;
			this->Success = false;
		}
	}
	if (!this->Success)
		LastError = serr.str();
	return this->Success;
}

void Configuration::printHelpMessage(stringstream &serr)
{
	serr << "[i] Coverage utility version " << VERSION << " (" << DATE << ")" << endl;
	serr << "[i] By " << AUTHOR << endl;
	serr << "[i] Usage: coverageUtil [arguments] <contigs.fasta> <read.coverage>" << endl;
	serr << "[i] -help                                               Print this message and exit." << endl;
	serr << "[i] -depth <filename>                                   Output per-base coverage depth in file <filename>. [off]" << endl;
}
