/*
 * dataLinker : creates abstract contig links from the available information sources.
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
	InputFileName = "";
	OutputFileName = "output.opt";
        ReadCoverageFileName = "";
	MaximumLinkHits = 5;
	NoOverlapDeviation = 0;
}

// Parses command line arguments. Returns true if successful.
bool Configuration::ProcessCommandLine(int argc, char *argv[])
{
	this->Success = true;
	stringstream serr;

	double weight = 1;
	int mapQ = 0;
	int minReadLength = 0;
	int maxEditDistance = 10000;

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
			else if (!strcmp("-mapq", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -mapq: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool newMapQSuccess;
				int newMapQ = Helpers::ParseInt(argv[i], newMapQSuccess);
				if (!newMapQSuccess)
				{
					serr << "[-] Parsing error in -mapq: weight must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				mapQ = newMapQ;
			}
			else if (!strcmp("-minlength", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -minlength: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool newMinReadLengthSuccess;
				int newMinReadlength = Helpers::ParseInt(argv[i], newMinReadLengthSuccess);
				if (!newMinReadLengthSuccess)
				{
					serr << "[-] Parsing error in -minlength: weight must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				minReadLength = newMinReadlength;
			}
			else if (!strcmp("-maxedit", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -maxedit: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				bool newMaxEditDistanceSuccess;
				int newMaxEditDistance = Helpers::ParseInt(argv[i], newMaxEditDistanceSuccess);
				if (!newMaxEditDistanceSuccess)
				{
					serr << "[-] Parsing error in -maxedit: weight must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				maxEditDistance = newMaxEditDistance;
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
			else if (!strcmp("-nooverlapdeviation", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -nooverlapdeviation: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				double newNoOverlapDeviation = atof(argv[i]);
				if (newNoOverlapDeviation < Helpers::Eps)
				{
					serr << "[-] Parsing error in -nooverlapdeviation: deviation must be a positive number." << endl;
					this->Success = false;
					break;
				}
				NoOverlapDeviation = newNoOverlapDeviation;
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
				this->PairedReadInputs.push_back(PairedInput(leftFileName, rightFileName, mu, sigma, false, weight, mapQ, minReadLength, maxEditDistance));
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
				this->PairedReadInputs.push_back(PairedInput(leftFileName, rightFileName, mu, sigma, true, weight, mapQ, minReadLength, maxEditDistance));
			}
                        else if (!strcmp("-seq", argv[i]))
			{
				if (argc - i - 1 < 2)
				{
					serr << "[-] Parsing error in -seq: must have 2 arguments." << endl;
					this->Success = false;
					break;
				}
				i++;
				string fileName = argv[i]; i++;
				bool sigmaSuccess;
				int sigma = Helpers::ParseInt(argv[i], sigmaSuccess);
				if (!sigmaSuccess || sigma < 0)
				{
					serr << "[-] Parsing error in -seq: sigma must be a non-negative number." << endl;
					this->Success = false;
					break;
				}
				this->SequenceInputs.push_back(SequenceInput(fileName, sigma, weight, minReadLength));
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
                        else if (!strcmp("-readcoverage", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -readcoverage: must have an argument." << endl;
					this->Success = false;
					break;
				}
				i++;
				this->ReadCoverageFileName = argv[i];
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
				this->BWAConfig.TmpPath = this->NovoAlignConfig.TmpPath = this->SAMToolsConfig.TmpPath = this->MummerTilerConfig.TmpPath = argv[i];
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
        serr << "[i] -minlength <length>                                 Set minimum length cutoff for information sources coming after the switch. [0]" << endl;
        serr << endl;
	serr << "[i] -mapq <score>                                       Set MapQ cutoff for information sources coming after the switch. [0]" << endl;
	serr << "[i] -maxedit <distance>                                 Set maximum edit distance cutoff for information sources coming after the switch. [10000]" << endl;
	serr << "[i] -maxhits <num>                                      Maximum number of allowed link hits. If a link has more hits, it is disregarded. [5]" << endl;
	serr << "[i] -nooverlapdeviation <num>                           Maximum allowed deviation from mean insert size when no overlaps are allowed. [disabled]" << endl;
	serr << "[i] -454 <left.fq> <right.fq> <mu> <sigma>              Process 454 paired reads with insert size <mu>+/-<sigma> into linking information." << endl;
	serr << "[i] -illumina <left.fq> <right.fq> <mu> <sigma>         Process Illumina paired reads with insert size <mu>+/-<sigma> into linking information." << endl;
        serr << endl;
        serr << "[i] -seq <reference.fa> <sigma>                         Process related sequences into linking information with <sigma> as standard deviation." << endl;
        serr << endl;
        serr << "[i] -readcoverage <filename>                            Produce contig read coverage data and output it to file <filename>. [disabled]" << endl;
	serr << "[i] -output <filename>                                  Output filename for optimzation information. [output.opt]" << endl;
	serr << "[i] -tmp <path>                                         Define scrap path for temporary files. [/tmp]" << endl;
	serr << "[i] BWA configuration options:" << endl;
	serr << "[i] -bwathreads <n>                                     Number of threads used in BWA alignment. [8]" << endl;
	serr << "[i] -bwahits <n>                                        Maximum number of alignment hits BWA should report. [1000]" << endl;
	serr << "[i] -bwaexact <yes/no>                                  Use exact matching in BWA? [no]";
        serr << endl;
}
