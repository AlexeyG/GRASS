/*
 * scaffoldOptimizer : solves the MIQP optimization and produces linear scaffold
 * sequences.
 * Copyright (C) 2011  Alexey Gritsenko
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 * 
 * 
 * Email: a.gritsenko@tudelft.nl
 * Mail: Delft University of Technology
 *       Faculty of Electrical Engineering, Mathematics, and Computer Science
 *       Department of Mediamatics
 *       P.O. Box 5031
 *       2600 GA, Delft, The Netherlands
 */

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
	RemoveAmbiguous = true;
	Sort = true;
	Bundle = true;
	BundlePerGroup = false;
	BundleAmbiguous = true;
	PrintMatrix = false;
	BundleDistance = 3.0;
	Erosion = 5.0;
        ExpectedCoverage = 0.0;
        UniquenessFCutoff = 5.0;
	InputFileName = "";
        ReadCoverageFileName = "";
	OutputFileName = "scaffold.fasta";
	SolutionOutputFileName = "";
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
			else if (!strcmp("-ambiguous", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -ambiguous: must have an argument." << endl;
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
					cerr << "[-] Parsing error in -ambiguous: argument must be yes/no." << endl;
					this->Success = false;
					break;
				}
				RemoveAmbiguous = !sw;
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
                        else if (!strcmp("-repeat-coverage", argv[i]))
			{
				if (argc - i - 1 < 2)
				{
					cerr << "[-] Parsing error in -repeat-coverage: must have at least two arguments argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool expectedCoverageSuccess;
				int expectedCoverage = Helpers::ParseInt(argv[i], expectedCoverageSuccess);
				if (!expectedCoverageSuccess || expectedCoverage <= 0)
				{
					cerr << "[-] Parsing error in -repeat-coverage: expected coverage must be a positive number." << endl;
					this->Success = false;
					break;
				}
				ExpectedCoverage = expectedCoverage;
                                i++;
                                ReadCoverageFileName = argv[i];
                                if (argc - i - 1 >= 1 && Helpers::IsNumber(argv[i + 1]))
                                {
                                    i++;
                                    bool uniquenessFCutoffSuccess;
                                    int uniquenessFCutoff = Helpers::ParseInt(argv[i], uniquenessFCutoffSuccess);
                                    if (!uniquenessFCutoffSuccess || uniquenessFCutoff <= 0)
                                    {
                                        cerr << "[i] Parsing error in -repeat-coverage: uniqueness F cutoff must be a position number." << endl;
                                        this->Success = false;
                                        break;
                                    }
                                    UniquenessFCutoff = uniquenessFCutoff;
                                }
			}
			else if (!strcmp("-erosion", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -erosion: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool erosionSuccess;
				int erosion = Helpers::ParseInt(argv[i], erosionSuccess);
				if (!erosionSuccess || erosion < 0)
				{
					cerr << "[-] Parsing error in -erosion: distance must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				Erosion = erosion;
			}
                        else if (!strcmp("-overlap-deviation", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -overlap-deviation: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool deviationSuccess;
				int deviation = Helpers::ParseInt(argv[i], deviationSuccess);
				if (!deviationSuccess || deviation < 0)
				{
					cerr << "[-] Parsing error in -overlap-deviation: length must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				OverlapperOptions.InitialOverlapDeviation = deviation;
			}
                        else if (!strcmp("-max-alignment", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -max-alignment: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool maxAlignemntSuccess;
				int maxAlignemnt = Helpers::ParseInt(argv[i], maxAlignemntSuccess);
				if (!maxAlignemntSuccess || maxAlignemnt < 0)
				{
					cerr << "[-] Parsing error in -max-alignment: length must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				OverlapperOptions.MaximumAlignmentLength = maxAlignemnt;
			}
                        else if (!strcmp("-no-split", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -no-split: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool noSplitSuccess;
				int noSplit = Helpers::ParseInt(argv[i], noSplitSuccess);
				if (!noSplitSuccess || noSplit < 0)
				{
					cerr << "[-] Parsing error in -no-split: length must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				OverlapperOptions.NoSplitOverlapLength = noSplit;
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
			else if (!strcmp("-solution-output", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -solution-output: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				SolutionOutputFileName = argv[i];
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
			else if (!strcmp("-ga-limit", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -ga-limit: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool timeLimitSuccess;
				int timeLimit = Helpers::ParseInt(argv[i], timeLimitSuccess);
				if (!timeLimitSuccess || timeLimit < 0)
				{
					cerr << "[-] Parsing error in -ga-limit: time must be a positive number." << endl;
					this->Success = false;
					break;
				}
				Options.GATimeLimit = timeLimit;
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
				int sw = 0;
				if (!strcasecmp(argv[i], "yes"))
					sw = 1;
				else if (!strcasecmp(argv[i], "no"))
					sw = 0;
				else if (!strcasecmp(argv[i], "more"))
					sw = 2;
				else
				{
					cerr << "[-] Parsing error in -verbose: argument must be yes/no/more." << endl;
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
	serr << "[i] -ambiguous <yes/no>                                 Use ambiguous contig links in optimization formulation. [no]" << endl;
	serr << "[i] -sort <yes/no>                                      Sort contig links in order to reduce non-zero elements in the optimization matrix. [yes]" << endl;
	serr << "[i] -bundle <yes/no>                                    Bundle contig links together? [yes]" << endl;
	serr << "[i] -bundle-groups <yes/no>                             Bundle contig links from different link groups? [yes]" << endl;
	serr << "[i] -bundle-ambiguous <yes/no>                          Bundle ambiguous and non-ambiguous contig links together? [yes]" << endl;
	serr << "[i] -bundle-distance <distance>                         Bundle contig links withing <distance> standard deviation from the median. [3]" << endl;
        serr << "[i] -repeat-coverage <exp. cov.> <cov. filename> [F]    Detect repeats using expected coverage and read coverage provided in a file. [5]" << endl;
	serr << "[i] -erosion <weight>                                   Remove contig links with weight smaller than <weigth> (should be used only with link bundling). [5]" << endl;
        serr << endl;
        serr << "[i] -overlap-deviation <length>                         Length by which contig distance is reduced to find better overlaps. [100]" << endl;
        serr << "[i] -max-alignment <length>                             Maximum contig overlap to perform global alignment on. Overlap deviation is not taken into account. [1500]" << endl;
        serr << "[i] -no-split <length>                                  Maximum predicted overlap length that is not confirmed by alignment, which does not cause scaffold splitting. [50]" << endl;
        serr << endl;
	serr << "[i] -time-limit <seconds>                               Time limit for a single run of CPLEX or GA in seconds. [infinite]" << endl;
	serr << "[i] -threads <n>                                        Number of threads a solver can use. [automatic]" << endl;
	serr << "[i] -cplex-opportunistic <yes/no>                       Use CPLEX opportunistic optimization mode. [yes]" << endl;
	serr << "[i] -cplex-heuristic <yes/no>                           Use CPLEX objective function heuristic. [yes]" << endl;
	serr << "[i] -cplex-suppress <yes/no>                            Suppress CPLEX output. [yes]" << endl;
	serr << "[i] -lp-limit <seconds>                                 Time in secconds for solving a single fixed optimization problem. [30]" << endl;
	serr << "[i] -lp-threads <seconds>                               Number of threads used for solving a single fixed optimization problem. [automatic]" << endl;
	serr << "[i] -lp-attempts <number>                               Number of attempts to solve a single fixed optimization problem. [3]" << endl;
	serr << "[i] -ga-limit <seconds>                                 Time in seconds for solving a single GA optimization problem. [unlimited]" << endl;
	serr << "[i] -ga-restarts <number>                               Number of restarts before exiting GA optimization. [4]" << endl;
	serr << "[i] -verbose <yes/no/more>                              Verbose output of solvers? [no]" << endl;
	serr << "[i] -output <output filename>                           Output filename for final scaffolds. [scaffold.fasta]" << endl;
	serr << "[i] -solution-output <output filename>                  Output filename for optimzation solution. [not output]" << endl;
}
