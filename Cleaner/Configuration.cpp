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
			else if (!strcmp("-454", argv[i]))
			{
				if (argc - i - 1 < 3)
				{
					serr << "[-] Parsing error in -454: must have 3 arguments." << endl;
					this->Success = false;
					break;
				}
				i++;
				string leftFileName = argv[i];
				i++;
				string rightFileName = argv[i];
				i++;
				string outputPrefix = argv[i];
				this->PairedReadInputs.push_back(PairedInput(leftFileName, rightFileName, outputPrefix, false));
			}
			else if (!strcmp("-illumina", argv[i]))
			{
				if (argc - i - 1 < 3)
				{
					serr << "[-] Parsing error in -illumina: must have 3 arguments." << endl;
					this->Success = false;
					break;
				}
				i++;
				string leftFileName = argv[i];
				i++;
				string rightFileName = argv[i];
				i++;
				string outputPrefix = argv[i];
				this->PairedReadInputs.push_back(PairedInput(leftFileName, rightFileName, outputPrefix, true));
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
				this->BWAConfig.TmpPath = this->NovoAlignConfig.TmpPath = this->SAMToolsConfig.TmpPath = argv[i];
			}
			else if (!strcmp("-bwathreads", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -bwathreads: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool threadsSuccess;
				BWAConfig.NumberOfThreads = Helpers::ParseInt(argv[i], threadsSuccess);
				if (!threadsSuccess || BWAConfig.NumberOfThreads <= 0)
				{
					serr << "[-] Parsing error in -bwathreads: number of threads must be a positive number." << endl;
					this->Success = false;
					break;
				}
			}
			else if (!strcmp("-bwahits", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -bwahits: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool hitsSuccess;
				BWAConfig.MaximumHits = Helpers::ParseInt(argv[i], hitsSuccess);
				if (!hitsSuccess || BWAConfig.MaximumHits <= 0)
				{
					serr << "[-] Parsing error in -bwahits: number of hits must be a positive number." << endl;
					this->Success = false;
					break;
				}
			}
			else if (!strcmp("-bwaexact", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -bwaexact: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool sw = false;
				if (!strcasecmp(argv[i], "yes"))
					sw = true;
				else if (!strcasecmp(argv[i], "no"))
					sw = false;
				else
				{
					cerr << "[-] Parsing error in -bwaexact: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				BWAConfig.ExactMatch = sw;
			}
			else if (i == argc - 1)
				this->ReferenceFileName = argv[argc - 1];
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
			serr << "[-] No reference file specified." << endl;
			this->Success = false;
		}
	}
	if (!this->Success)
		LastError = serr.str();
	return this->Success;
}

void Configuration::printHelpMessage(stringstream &serr)
{
	serr << "[i] Repetitive read filter " << VERSION << " (" << DATE << ")" << endl;
	serr << "[i] By " << AUTHOR << endl;
	serr << "[i] Usage: dataLinker [arguments] <referece.fasta>" << endl;
	serr << "[i] -help                                               Print this message and exit." << endl;
	serr << "[i] -454 <left.fq> <right.fq> <output prefix>           Process 454 paired reads and output the filtered reads with new prefix." << endl;
	serr << "[i] -illumina <left.fq> <right.fq> <output prefix>      Process Illumina paired reads and output the filtered reads with new prefix." << endl;
	serr << "[i] -tmp <path>                                         Define scrap path for temporary files. [/tmp]" << endl;
	serr << "[i] BWA configuration options:" << endl;
	serr << "[i] -bwathreads <n>                                     Number of threads used in BWA alignment. [8]" << endl;
	serr << "[i] -bwahits <n>                                        Maximum number of alignment hits BWA should report. [1000]" << endl;
	serr << "[i] -bwaexact <yes/no>                                  Use exact matching in BWA? [no]";
}
