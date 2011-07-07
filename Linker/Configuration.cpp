#include "Configuration.h"
#include "Globals.h"
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
	InputFileName = "";
	OutputFileName = "output.opt";
	MaximumLinkHits = 5;
}

// Parses command line arguments. Returns true if successful.
bool Configuration::ProcessCommandLine(int argc, char *argv[])
{
	this->Success = true;
	stringstream serr;

	double weight = 1;

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
			else if (!strcmp("-weight", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -weight: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				double newWeight = atof(argv[i]);
				if (newWeight < Helpers::Eps)
				{
					serr << "[-] Parsing error in -weight: weight must be a positive number." << endl;
					this->Success = false;
					break;
				}
				weight = newWeight;
			}
			else if (!strcmp("-maxhits", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -maxhits: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool hitsSuccess;
				MaximumLinkHits = Helpers::ParseInt(argv[i], hitsSuccess);
				if (!hitsSuccess || MaximumLinkHits <= 0)
				{
					serr << "[-] Parsing error in -maxhits: number of hits must be a positive number." << endl;
					this->Success = false;
					break;
				}
			}
			else if (!strcmp("-454", argv[i]))
			{
				if (argc - i - 1 < 4)
				{
					serr << "[-] Parsing error in -454: must have 4 arguments." << endl;
					this->Success = false;
					break;
				}
				i++;
				string leftFileName = argv[i]; i++;
				string rightFileName = argv[i]; i++;
				bool muSuccess;
				int mu = Helpers::ParseInt(argv[i], muSuccess);
				if (!muSuccess || mu < 0)
				{
					serr << "[-] Parsing error in -454: mu must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool sigmaSuccess;
				int sigma = Helpers::ParseInt(argv[i], sigmaSuccess);
				if (!sigmaSuccess || sigma < 0)
				{
					serr << "[-] Parsing error in -454: sigma must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				this->PairedReadInputs.push_back(PairedInput(leftFileName, rightFileName, mu, sigma, false, weight));
			}
			else if (!strcmp("-illumina", argv[i]))
			{
				if (argc - i - 1 < 4)
				{
					serr << "[-] Parsing error in -illumina: must have 4 arguments." << endl;
					this->Success = false;
					break;
				}
				i++;
				string leftFileName = argv[i]; i++;
				string rightFileName = argv[i]; i++;
				bool muSuccess;
				int mu = Helpers::ParseInt(argv[i], muSuccess);
				if (!muSuccess || mu < 0)
				{
					serr << "[-] Parsing error in -illumina: mu must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool sigmaSuccess;
				int sigma = Helpers::ParseInt(argv[i], sigmaSuccess);
				if (!sigmaSuccess || sigma < 0)
				{
					serr << "[-] Parsing error in -illumina: sigma must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				this->PairedReadInputs.push_back(PairedInput(leftFileName, rightFileName, mu, sigma, true, weight));
			}
			else if (!strcmp("-output", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -output: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				this->OutputFileName = argv[i];
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
	serr << "[i] Contig linker version " << VERSION << " (" << DATE << ")" << endl;
	serr << "[i] By " << AUTHOR << endl;
	serr << "[i] Usage: dataLinker [arguments] <sequence.fasta>" << endl;
	serr << "[i] -help                                               Print this message and exit." << endl;
	serr << "[i] -weight <weight>                                    Set weight for information sources coming after the switch. [1]" << endl;
	serr << "[i] -maxhits <num>                                      Maximum number of allowed link hits. If a link has more hits, it is disregarded. [5]" << endl;
	serr << "[i] -454 <left.fq> <right.fq> <mu> <sigma>              Process 454 paired reads with insert size <mu>+/-<sigma> into linking information." << endl;
	serr << "[i] -illumina <left.fq> <right.fq> <mu> <sigma>         Process Illumina paired reads with insert size <mu>+/-<sigma> into linking information." << endl;
	serr << "[i] -output [output filename]                           Output filename for optimzation information. [output.opt]" << endl;
	serr << "[i] -tmp <path>                                         Define scrap path for temporary files. [/tmp]" << endl;
	serr << "[i] BWA configuration options:" << endl;
	serr << "[i] -bwathreads <n>                                     Number of threads used in BWA alignment. [8]" << endl;
	serr << "[i] -bwahits <n>                                        Maximum number of alignment hits BWA should report. [1000]" << endl;
}

// Construtor with default configuration parameter settings.
BWAConfiguration::BWAConfiguration()
{
	TmpPath = "/tmp";
	NumberOfThreads = 8;
	MaximumHits = 1000;
	IndexCommand = "bwa index -a is -p %s %s >& /dev/null";
	SuffixArrayCommand = "bwa aln -t %i -f %s %s %s >& /dev/null";
	AlignSingleEndCommand = "bwa samse -f %s -n %i %s %s %s >& /dev/null";
}

// Construtor with default configuration parameter settings.
NovoAlignConfiguration::NovoAlignConfiguration()
{
	TmpPath = "/tmp";
	IndexCommand = "novoindex -m %s %s >& /dev/null";
	AlignSingleEndCommand = "(novoalign -d %s -f %s -o SAM -r All > %s) >& /dev/null";
}

// Construtor with default configuration parameter settings.
SAMToolsConfiguration::SAMToolsConfiguration()
{
	TmpPath = "/tmp";
	ConvertCommand = "samtools view -b -S -o %s %s >& /dev/null";
}
