/*
 * dataSelector : a support tool used in debugging to create simulated scaffolding
 * problems given a reference sequence and paired reads.
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
#include "Helpers.h"
#include "Defines.h"
#include <sstream>
#include <cstring>
#include <cstdlib>

// Configuration segment class ctor
Segment::Segment(int start, int finish, int chromosome) : Start(start), Finish(finish), Chromosome(chromosome)
{
}

// Comparison operator
bool Segment::operator< (const Segment &other) const
{
	if (Chromosome < other.Chromosome)
		 return true;
	if (Chromosome > other.Chromosome)
		return false;
	if (Start < other.Start)
		return true;
	if (Start > other.Start)
		return false;
	if (Finish > other.Finish)
		return true;
	if (Finish < other.Finish)
		return false;
	return false;
}

// Configuration selection class ctor
ConfigSelect::ConfigSelect(int length, int chromosome, int count) : Length(length), Chromosome(chromosome), Count(count)
{
}

// Configuration input bam class ctor
PairedBam::PairedBam(const string &inputBam1, const string &inputBam2, const string &prefix) : InputBam1(inputBam1), InputBam2(inputBam2), OutputPrefix(prefix)
{
}

PairedSimulation::PairedSimulation(double readLengthMean, double readLengthStd, double insertSizeMean, double insertSizeStd, double depth, bool isIllumina, const string &prefix)
	: ReadLengthMean(readLengthMean), ReadLengthStd(readLengthStd), InsertSizeMean(insertSizeMean), InsertSizeStd(insertSizeStd), Depth(depth), IsIllumina(isIllumina), OutputPrefix(prefix)
{
}

// Configuration class defualt ctor
Configuration::Configuration()
{
	InputFastaFileName = "";
	OutputFastaFileName = "";
	PrintChromosomeInfo = false;
	Select.clear();
	Segments.clear();
	PairedAlignment.clear();
}

bool Configuration::ProcessCommandLine(int argc, char *argv[])
{
	Success = true;
	stringstream serr;
	if (argc == 1)
	{
		serr << "[-] Not enough arguments. Consult -help." << endl;
		Success = false;
	}
	else
	{
		OutputFastaFileName = "out.fasta";
		int i = 1;
		while (i < argc)
		{
			if (!strcmp("-help", argv[i]) || !strcmp("-h", argv[i]))
			{
				printHelpMessage(serr);
				Success = false;
				break;
			}
			else if (!strcmp("-chromosomes", argv[i]))
				PrintChromosomeInfo = true;
			else if (!strcmp("-segment", argv[i]))
			{
				if (argc - i - 1 < 2)
				{
					serr << "[-] Parsing error in -segment: must have at least two arguments." << endl;
					Success = false;
					break;
				}
				i++;
				int start = atoi(argv[i]);
				if (start <= 0)
				{
					serr << "[-] Parsing error in -segment: start must be a positive number." << endl;
					Success = false;
					break;
				}
				i++;
				int length = atoi(argv[i]); 
				if (length <= 0)
				{
					serr << "[-] Parsing error in -segment: length must be a positive number." << endl;
					Success = false;
					break;
				}
				int chromosome = 1;
				if (i + 1 < argc && Helpers::IsNumber(argv[i + 1]))
				{
					i++;
					chromosome = atoi(argv[i]);
				}
				if (chromosome <= 0)
				{
					serr << "[-] Parsing error in -segment: chromosome must be a positive number." << endl;
					Success = false;
					break;
				}
				start--;
				Segments.push_back(Segment(start, start + length - 1, chromosome));
			}
			else if (!strcmp("-select", argv[i]))
			{
				if (argc - i - 1 < 2)
				{
					serr << "[-] Parsing error in -select: must have at least two argument." << endl;
					Success = false;
					break;
				}
				i++;
				int selectCount = 1;
				int selectChromosome = atoi(argv[i + 1]);
				int selectLength = atoi(argv[i]);
				if (selectLength <= 0)
				{
					serr << "[-] Parsing error in -select: length must be a positive number." << endl;
					Success = false;
					break;
				}
				if (selectChromosome <= 0)
				{
					serr << "[-] Parsing error in -select: chromosome # must be a positive number." << endl;
					Success = false;
					break;
				}
				i++;
				if (i + 1 < argc && Helpers::IsNumber(argv[i + 1]))
				{
					i++;
					selectCount = atoi(argv[i]);
				}
				if (selectCount <= 0)
				{
					serr << "[-] Parsing error in -select: count must be a positive number." << endl;
					Success = false;
					break;
				}
				Select.push_back(ConfigSelect(selectLength, selectChromosome, selectCount));
			}
			else if (!strcmp("-selectall", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -selectall: must have at least one argument." << endl;
					Success = false;
					break;
				}
				i++;
				int selectCount = 1;
				int selectChromosome = -1;
				int selectLength = atoi(argv[i]);
				if (selectLength <= 0)
				{
					serr << "[-] Parsing error in -selectall: length must be a positive number." << endl;
					Success = false;
					break;
				}
				if (i + 1 < argc && Helpers::IsNumber(argv[i + 1]))
				{
					i++;
					selectCount = atoi(argv[i]);
				}
				if (selectCount <= 0)
				{
					serr << "[-] Parsing error in -selectall: count must be a positive number." << endl;
					Success = false;
					break;
				}
				Select.push_back(ConfigSelect(selectLength, selectChromosome, selectCount));
			}
			else if (!strcmp("-output", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					serr << "[-] Parsing error in -output: must have an argument." << endl;
					Success = false;
					break;
				}
				i++;
				OutputFastaFileName = argv[i];
			}
			else if (!strcmp("-paired", argv[i]))
			{
				if (argc - i - 1 < 3)
				{
					serr << "[-] Parsing error in -paired: must have three arguments." << endl;
					Success = false;
					break;
				}
				i++;
				string inputBam1 = argv[i];
				i++;
				string inputBam2 = argv[i];
				i++;
				string outputPrefix = argv[i];
				PairedAlignment.push_back(PairedBam(inputBam1, inputBam2, outputPrefix));
			}
			else if (!strcmp("-simulate-paired", argv[i]))
			{
				if (argc - i - 1 < 5)
				{
					serr << "[-] Parsing error in -simulate-paired: must have at least five arguments." << endl;
					Success = false;
					break;
				}
				bool readLengthMeanSuccess;
				i++; int readLengthMean = Helpers::ParseInt(argv[i], readLengthMeanSuccess);
				if (readLengthMean < 0 || !readLengthMeanSuccess)
				{
					serr << "[-] Parsing error in -simulate-paired: mean read length must be a positive number." << endl;
					Success = false;
					break;
				}
				bool readLengthStdSuccess;
				i++; int readLengthStd = Helpers::ParseInt(argv[i], readLengthStdSuccess);
				if (readLengthStd < 0 || !readLengthStdSuccess)
				{
					serr << "[-] Parsing error in -simulate-paired: read length deviation must be a positive number." << endl;
					Success = false;
					break;
				}
				bool insertMeanSuccess;
				i++; int insertMean = Helpers::ParseInt(argv[i], insertMeanSuccess);
				if (insertMean < 0 || !insertMeanSuccess)
				{
					serr << "[-] Parsing error in -simulate-paired: mean insert size must be a positive number." << endl;
					Success = false;
					break;
				}
				bool insertStdSuccess;
				i++; int insertStd = Helpers::ParseInt(argv[i], insertStdSuccess);
				if (insertStd < 0 || !insertStdSuccess)
				{
					serr << "[-] Parsing error in -simulate-paired: insert size deviation must be a positive number." << endl;
					Success = false;
					break;
				}
				bool depthSuccess;
				i++; int depth = Helpers::ParseInt(argv[i], depthSuccess);
				if (depth < 0 || !depthSuccess)
				{
					serr << "[-] Parsing error in -simulate-paired: depth must be a positive number." << endl;
					Success = false;
					break;
				}
				i++; string prefix = argv[i];
				bool isIllumina = true;
				if (i + 1 < argc && (strcasecmp(argv[i + 1], "illumina") == 0 || strcasecmp(argv[i + 1], "454") == 0))
				{
					i++;
					isIllumina = strcasecmp(argv[i], "illumina") == 0;
				}
				PairedReadSimulation.push_back(PairedSimulation(readLengthMean, readLengthStd, insertMean, insertStd, depth, isIllumina, prefix));
			}
			else if (i == argc - 1)
				InputFastaFileName = argv[argc - 1];
			else
			{
				serr << "[i] Unknown argument: " << argv[i] << endl;
				Success = false;
				break;
			}
			i++;
		}

		if (InputFastaFileName == "")
		{
			serr << "[-] No sequence file given." << endl;
			Success = false;
		}
	}

	if (!Success)
		LastError = serr.str();
	return Success;
}

void Configuration::printHelpMessage(stringstream &serr)
{
	serr << "[i] Data selector version " << VERSION << " (" << DATE << ")" << endl;
	serr << "[i] By " << AUTHOR << endl;
	serr << "[i] Usage: dataSelector [arguments] <sequence.fasta>" << endl;
	serr << "[i] -help                                                     Print this message and exit." << endl;
	serr << "[i] -chromosomes                                              Print chromosome numbers and exit." << endl;
	serr << "[i] -segment <start> <length> [chromosome #]                  Selects segment [start; start + length) from given chromosome. [1]" << endl;
	serr << "                                                              Indices are 1-based." << endl;
	serr << "[i] -select <length> <chromosome #> [count]                   Randomly selects [count] regions of given length in the given chromosome." << endl;
	serr << "[i] -selectall <length> [count]                               Randomly selects [count] regions of given length in every chromosome." << endl;
	serr << "[i] -paired <input1.bam> <input2.bam> <output prefix>         Processes paired end BAM alignments and outputs two read files with given prefix." << endl;
	serr << "[i] -simulate-paired <mean-length> <std-length>               Simulate paired reads with given read length and insert size and output two read" << endl;
	serr << "    <mean-insert> <std-insert> <depth> <output prefix> [type] files with given prefix. Type is Illumina or 454. [Illumina]" << endl;
	serr << "[i] -output [sequence.fasta]                                  Output filename for the selected sequences. [out.fasta]" << endl;
}
