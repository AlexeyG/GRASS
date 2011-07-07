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
	Sort = true;
	Bundle = true;
	BundlePerGroup = false;
	BundleAmbiguous = true;
	PrintMatrix = false;
	BundleDistance = 3.0;
	InputFileName = "";
	OutputFileName = "scaffold.fasta";
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
			else if (!strcmp("-sort", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -sort: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -sort: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				Sort = sw;
			}
			else if (!strcmp("-bundle", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -bundle: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -bundle: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				Bundle = sw;
			}
			else if (!strcmp("-bundle-groups", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -bundle-groups: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -bundle-groups: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				BundlePerGroup = !sw;
			}
			else if (!strcmp("-bundle-ambiguous", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -bundle-ambiguous: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -bundle-ambiguous: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				BundleAmbiguous = sw;
			}
			else if (!strcmp("-bundle-distance", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -bundle-distance: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool distanceSuccess;
				int distance = Helpers::ParseInt(argv[i], distanceSuccess);
				if (!distanceSuccess || distance <= 0)
				{
					cerr << "[-] Parsing error in -bundle-distance: distance must be a positive number." << endl;
					this->Success = false;
					break;
				}
				BundleDistance = distance;
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
				OutputFileName = argv[i];
			}
			else if (!strcmp("-print-matrix", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -print-matrix: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -print-matrix: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				PrintMatrix = sw;
			}
			else if (!strcmp("-time-limit", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -time-limit: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool timeLimitSuccess;
				int timeLimit = Helpers::ParseInt(argv[i], timeLimitSuccess);
				if (!timeLimitSuccess || timeLimit < 0)
				{
					cerr << "[-] Parsing error in -time-limit: time must be a positive number." << endl;
					this->Success = false;
					break;
				}
				Options.TimeLimit = timeLimit;
			}
			else if (!strcmp("-threads", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -threads: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool threadsSuccess;
				int threads = Helpers::ParseInt(argv[i], threadsSuccess);
				if (!threadsSuccess || threads < 0)
				{
					cerr << "[-] Parsing error in -threads: number of threads must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				Options.Threads = threads;
			}
			else if (!strcmp("-cplex-opportunistic", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -cplex-opportunistic: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -cplex-opportunistic: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				Options.UseOpportunisticSearch = sw;
			}
			else if (!strcmp("-cplex-heuristic", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -cplex-heuristic: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -cplex-heuristic: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				Options.UseObjectiveHeuristic = sw;
			}
			else if (!strcmp("-cplex-suppress", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -cplex-suppress: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -cplex-suppress: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				Options.SuppressOutput = sw;
			}
			else if (!strcmp("-lp-limit", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -lp-limit: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool timeLimitSuccess;
				int timeLimit = Helpers::ParseInt(argv[i], timeLimitSuccess);
				if (!timeLimitSuccess || timeLimit < 0)
				{
					cerr << "[-] Parsing error in -lp-limit: time must be a positive number." << endl;
					this->Success = false;
					break;
				}
				Options.LPTimeLimit = timeLimit;
			}
			else if (!strcmp("-lp-threads", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -lp-threads: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool threadsSuccess;
				int threads = Helpers::ParseInt(argv[i], threadsSuccess);
				if (!threadsSuccess || threads < 0)
				{
					cerr << "[-] Parsing error in -lp-threads: number of threads must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				Options.LPThreads = threads;
			}
			else if (!strcmp("-lp-attempts", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -lp-attempts: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool attemptsSuccess;
				int attempts = Helpers::ParseInt(argv[i], attemptsSuccess);
				if (!attemptsSuccess || attempts < 0)
				{
					cerr << "[-] Parsing error in -lp-attempts: number of attempts must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				Options.LPAttempts = attempts;
			}
			else if (!strcmp("-ga-restarts", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -ga-restarts: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool restartsSuccess;
				int restarts = Helpers::ParseInt(argv[i], restartsSuccess);
				if (!restartsSuccess || restarts < 0)
				{
					cerr << "[-] Parsing error in -ga-restarts: number of attempts must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				Options.GARestarts = restarts;
			}
			else if (!strcmp("-verbose", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -verbose: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -verbose: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				Options.VerboseOutput = sw;
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
	serr << "[i] -sort <yes/no>                                      Sort contig links in order to reduce non-zero elements in the optimization matrix. [yes]" << endl;
	serr << "[i] -bundle <yes/no>                                    Bundle contig links together? [yes]" << endl;
	serr << "[i] -bundle-groups <yes/no>                             Bundle contig links from different link groups? [yes]" << endl;
	serr << "[i] -bundle-ambiguous <yes/no>                          Bundle ambiguous and non-ambiguous contig links together? [yes]" << endl;
	serr << "[i] -bundle-distance <distance>                         Bundle contig links withing <distance> standard deviation from the median. [3]" << endl;
	serr << "[i] -time-limit <seconds>                               Time limit for a single run of CPLEX or GA in seconds. [infinite]" << endl;
	serr << "[i] -threads <n>                                        Number of threads a solver can use. [automatic]" << endl;
	serr << "[i] -cplex-opportunistic <yes/no>                       Use CPLEX opportunistic optimization mode. [yes]" << endl;
	serr << "[i] -cplex-heuristic <yes/no>                           Use CPLEX objective function heuristic. [yes]" << endl;
	serr << "[i] -cplex-suppress <yes/no>                            Suppress CPLEX output. [yes]" << endl;
	serr << "[i] -lp-limit <seconds>                                 Time in secconds for solving a single fixed optimization problem. [30]" << endl;
	serr << "[i] -lp-threads <seconds>                               Number of threads used for solving a single fixed optimization problem. [automatic]" << endl;
	serr << "[i] -lp-attempts <number>                               Number of attempts to solve a single fixed optimization problem. [3]" << endl;
	serr << "[i] -ga-restarts <number>                               Number of restarts before exiting GA optimization. [unlimited]" << endl;
	serr << "[i] -verbose <yes/no>                                   Verbose output of solvers? [no] [3]" << endl;
	serr << "[i] -output [output filename]                           Output filename for optimzation information. [scaffold.fasta]" << endl;
}
