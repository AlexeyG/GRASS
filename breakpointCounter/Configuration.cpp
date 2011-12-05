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
	ScaffoldFileName = "";
        ReferenceFileName = "";
        MinBases = 90;
}

// Parses command line arguments. Returns true if successful.
bool Configuration::ProcessCommandLine(int argc, char *argv[])
{
    this->Success = true;
    stringstream serr;

    if (argc < 2)
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
            else if (!strcmp("-minbases", argv[i]))
            {
                if (argc - i - 1 < 1)
                {
                    serr << "[-] Parsing error in -minbases: must have an argument." << endl;
                    this->Success = false;
                    break;
                }
                i++;
                bool minBasesSuccess;
                int minBases = Helpers::ParseInt(argv[i], minBasesSuccess);
                if (!minBasesSuccess || minBases < 0)
                {
                    cerr << "[-] Parsing error in -minbases: number of bases must be a non-negative number." << endl;
                    this->Success = false;
                    break;
                }
            }
            else if (!strcmp("-tmp", argv[i]))
            {
                if (argc - i - 1 < 1)
                {
                    serr << "[-] Parsing error in -tmp: must have an argument." << endl;
                    this->Success = false;
                    break;
                }
                i++;
                MummerConfig.TmpPath = argv[i];
            }
            else if (i == argc - 2)
                this->ReferenceFileName = argv[argc - 2];
            else if (i == argc - 1)
                this->ScaffoldFileName = argv[argc - 1];
            else
            {
                serr << "[-] Unknown argument: " << argv[i] << endl;
                this->Success = false;
                break;
            }
            i++;
        }
        if (this->ReferenceFileName == "")
        {
            serr << "[-] No reference input file specified." << endl;
            this->Success = false;
        }
        if (this->ScaffoldFileName == "")
        {
            serr << "[-] No scaffolds input file specified." << endl;
            this->Success = false;
        }
    }
    if (!this->Success)
            LastError = serr.str();
    return this->Success;
}

void Configuration::printHelpMessage(stringstream &serr)
{
    serr << "[i] Scaffold breakpoint counter version " << VERSION << " (" << DATE << ")" << endl;
    serr << "[i] By " << AUTHOR << endl;
    serr << "[i] Usage: breakpointCounter [arguments] <reference.fasta> <scaffolds.fasta>" << endl;
    serr << "[i] -help                                               Print this message and exit." << endl;
    serr << "[i] -tmp <path>                                         Define scrap path for temporary files. [/tmp]" << endl;
    serr << "[i] -minbases <number>                                  Minimum number of aligned bases to take into account. [90]" << endl;
}
