#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include "Globals.h"
#include "Configuration.h"
#include "Sequence.h"
#include "ContigInformation.h"
#include "Reader.h"
#include "Writer.h"
#include "Helpers.h"

using namespace std;

Configuration configuration;
vector<FastASequence> contigs;
vector<ContigInformation> infos;

void printHelpMessage(void)
{
	cerr << "[i] Data simulator version " << VERSION << " (" << DATE << ")" << endl;
	cerr << "[i] By " << AUTHOR << endl;
	cerr << "[i] Usage: dataSimulator [arguments] <sequence.fasta>" << endl;
	cerr << "[i] -help                                               Print this message and exit." << endl;
	cerr << "[i] -split <mu> <sigma> [# splits = 1]                  Introduce [# splits] gaps of mean size <mu> and variance <sigma> in the input contigs." << endl;
	cerr << "[i] -limit <len>                                        Do not split the contig if the contig becomes shorter than <len> after split." << endl;
	cerr << "[i] -overlap <yes/no>                                   Allow contig overlaps (negative gap lengths will be treated as overlaps) [no]." << endl;
	cerr << "[i] -flip <yes/no>                                      Randomly flip orientation of the created contigs [yes]." << endl;
	cerr << "[i] -shuffle <yes/no>                                   Randomly shuffle created contigs [yes]." << endl;
	cerr << "[i] -output [output filename]                           Output filename for the simulated data. [out.fasta]" << endl;
}

Configuration processCommandLine(int argc, char *argv[])
{
	Configuration config;
	config.Success = true;

	if (argc == 1)
	{
		cerr << "[-] Not enough arguments. Consult -help." << endl;
		config.Success = false;
	}
	else
	{
		int i = 1;
		while (i < argc)
		{
			if (!strcmp("-help", argv[i]) || !strcmp("-h", argv[i]))
			{
				printHelpMessage();
				config.Success = false;
				break;
			}
			else if (!strcmp("-split", argv[i]))
			{
				if (argc - i - 1 < 2)
				{
					cerr << "[-] Parsing error in -split: must have at least two argument." << endl;
					config.Success = false;
					break;
				}
				i++;
				bool muSuccess;
				int mu = Helpers::ParseInt(argv[i], muSuccess);
				if (!muSuccess || mu < 0)
				{
					cerr << "[-] Parsing error in -split: mu must be a non-negative number." << endl;
					config.Success = false;
					break;
				}
				i++;
				bool sigmaSuccess;
				int sigma = Helpers::ParseInt(argv[i], sigmaSuccess);
				if (!sigmaSuccess || sigma < 0)
				{
					cerr << "[-] Parsing error in -split: sigma must be a non-negative number." << endl;
					config.Success = false;
					break;
				}
				int count = 1;
				bool countSuccess;
				if (i + 1 < argc && Helpers::IsNumber(argv[i + 1]))
				{
					i++;
					count = Helpers::ParseInt(argv[i], countSuccess);
					if (!countSuccess || count <= 0)
					{
						cerr << "[-] Parsing error in -split: count must be a positive number." << endl;
						config.Success = false;
						break;
					}
				}
				config.Splits.push_back(Split(mu, sigma, count));
			}
			else if (!strcmp("-limit", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -limit: must have an argument." << endl;
					config.Success = false;
					break;
				}
				i++;
				bool limitSuccess;
				int limit = Helpers::ParseInt(argv[i], limitSuccess);
				if (!limitSuccess || limit <= 0)
				{
					cerr << "[-] Parsing error in -limit: length must be a positive number." << endl;
					config.Success = false;
					break;
				}
				config.Limit = limit;
			}
			else if (!strcmp("-overlap", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -overlap: must have an argument." << endl;
					config.Success = false;
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
					cerr << "[-] Parsing error in -overlap: argument must be yes/no." << endl;
					config.Success = false;
					break;
				}
				config.AllowContigOverlap = sw;
			}
			else if (!strcmp("-flip", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -flip: must have an argument." << endl;
					config.Success = false;
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
					cerr << "[-] Parsing error in -flip: argument must be yes/no." << endl;
					config.Success = false;
					break;
				}
				config.FlipOrientation = sw;
			}
			else if (!strcmp("-shuffle", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -shuffle: must have an argument." << endl;
					config.Success = false;
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
					cerr << "[-] Parsing error in -shuffle: argument must be yes/no." << endl;
					config.Success = false;
					break;
				}
				config.ShuffleContigs = sw;
			}
			else if (!strcmp("-output", argv[i]))
			{
				if (argc - i - 1 < 1)
				{
					cerr << "[-] Parsing error in -output: must have an argument." << endl;
					config.Success = false;
					break;
				}
				i++;
				config.OutputFileName = argv[i];
			}
			else if (i == argc - 1)
				config.InputFileName = argv[argc - 1];
			else
			{
				cerr << "[-] Unknown argument: " << argv[i] << endl;
				config.Success = false;
				break;
			}
			i++;
		}
		if (config.InputFileName == "")
		{
			cerr << "[-] No input file specified." << endl;
			config.Success = false;
		}
	}

	return config;
}

bool readContigs(const string &inputFile, vector<FastASequence> &contigs, vector<ContigInformation> &infos)
{
	FastAReader reader;
	if (!reader.Open(inputFile))
		return false;
	
	int read = reader.Read(contigs);
	infos.resize(read);
	reader.Close();
	return read != 0;
}

bool writeContigs(const string &outputFile, const vector<FastASequence> &contigs)
{
	FastAWriter writer;
	if (!writer.Open(outputFile))
		return false;

	bool success = writer.Write(contigs);
	writer.Close();
	return success;
}

bool canIntroduceContigGap(const FastASequence &contig, int gap, int &coord, int limit, bool canOverlap)
{
	int len = contig.Nucleotides.size();
	if (len - gap < 2 * limit)
		return false;
	coord = rand() % (len - gap - 2 * limit + 1) + limit;
	return true;
}

bool canSplitContig(const FastASequence &contig, const Split &split, int limit, bool allowOverlap, int &gap, int &coord)
{
	for (int att = 0; att < MaxAttempts; att++)
	{
		gap = Helpers::RandomNormal(split.Mean, split.Variance);
		if (gap < 0 && !allowOverlap)
			continue;
		if (canIntroduceContigGap(contig, gap, coord, limit, allowOverlap))
			return true;
	}
	return false;
}

bool canSplit(vector<FastASequence> &contigs, const Split &split, int limit, bool allowOverlap, int &pos, int &gap, int &coord)
{
	int initialPos = pos;
	while (!canSplitContig(contigs[pos], split, limit, allowOverlap, gap, coord))
	{
		pos = (pos + 1) % contigs.size();
		if (pos == initialPos)
				return false;
	}
	return true;
}

pair<FastASequence, FastASequence> cutContig(const FastASequence &contig, int gap, int coord)
{
	// [0; sel- 1] [sel; sel + gap - 1] [sel + gap; len - 1]
	string seqA, seqB;
	int len = contig.Nucleotides.length();
	seqA = contig.Nucleotides.substr(0, coord);
	seqB = contig.Nucleotides.substr(coord + gap, len - coord - gap);
	FastASequence contigA = FastASequence(seqA, contig.Comment);
	FastASequence contigB = FastASequence(seqB, contig.Comment);

	return pair<FastASequence, FastASequence>(contigA, contigB);
}

pair<ContigInformation, ContigInformation> cutContigInformation(const ContigInformation &info, int gap, int coord)
{
	ContigInformation a = info;
	ContigInformation b(a.Position + coord + gap, a.ReverseOrientation);
	return pair<ContigInformation, ContigInformation>(a, b);
}

void splitContigs(const Configuration &config, vector<FastASequence> &contigs, vector<ContigInformation> &infos)
{
	int nSplit = config.Splits.size();
	int pos = 0;
	for (int i = 0; i < nSplit; i++)
	{
		int gap, coord;
		Split split = config.Splits[i];
		int cnt = split.Count;
		cerr << "[i] Processing splits of mean " << split.Mean << " and variance " << split.Variance << " (" << cnt << " splits)." << endl;
		while (cnt > 0 && canSplit(contigs, split, config.Limit, config.AllowContigOverlap, pos, gap, coord))
		{
			pair<FastASequence, FastASequence> cut = cutContig(contigs[pos], gap, coord);
			pair<ContigInformation, ContigInformation> cutInfo = cutContigInformation(infos[pos], gap, coord);
			contigs[pos] = cut.first;
			infos[pos] = cutInfo.first;
			contigs.insert(contigs.begin() + pos + 1, cut.second);
			infos.insert(infos.begin() + pos + 1, cutInfo.second);
			cnt--;
			cerr << "   [+] Introduced a gap of size " << gap << " at position " << coord << " (" << cnt << " left)." << endl;
			pos = (pos + 2) % contigs.size();
		}
		if (cnt > 0)
			cerr << "   [-] Unable to introduce all splits (" << cnt << " left)." << endl;
	}
}

void flipOrientation(vector<FastASequence> &contigs, vector<ContigInformation> &infos)
{
	int n = contigs.size();
	for (int i = 0; i < n; i++)
		if (rand() % 2 == 0)
		{
			contigs[i].ReverseCompelement();
			infos[i].ReverseOrientation = !infos[i].ReverseOrientation;
		}
}

void shuffleContigs(vector<FastASequence> &contigs, vector<ContigInformation> &infos)
{
	int n = contigs.size();
	if (n == 0)
		return;
	for (int i = 0; i < n; i++)
	{
		int j = rand() % n;
		FastASequence tmp = contigs[i];
		contigs[i] = contigs[j];
		contigs[j] = tmp;
		ContigInformation tmpInfo = infos[i];
		infos[i] = infos[j];
		infos[j] = tmpInfo;
	}
}

void updateContigComments(vector<FastASequence> &contigs, const vector<ContigInformation> &infos)
{
	int n = contigs.size();
	for (int i = 0; i < n; i++)
		contigs[i].Comment = infos[i].FormatName() + contigs[i].Comment;
}

int main(int argc, char *argv[])
{
	srand((unsigned int)time(NULL));

	configuration = processCommandLine(argc, argv);
	if (configuration.Success)
	{
		if (!readContigs(configuration.InputFileName, contigs, infos))
		{
			cerr << "[-] Unable to read contigs from file (" << configuration.InputFileName << ")." << endl;
			return -1;
		}
		cerr << "[+] Read input contigs from file (" << configuration.InputFileName << ")." << endl;
		splitContigs(configuration, contigs, infos);
		if (configuration.FlipOrientation)
		{
			flipOrientation(contigs, infos);
			cerr << "[+] Flipped contig orientation." << endl;
		}
		if (configuration.ShuffleContigs)
		{
			shuffleContigs(contigs, infos);
			cerr << "[+] Shuffled contigs." << endl;
		}
		updateContigComments(contigs, infos);
		if (!writeContigs(configuration.OutputFileName, contigs))
		{
			cerr << "[-] Unable to write contigs to file (" << configuration.OutputFileName << ")." << endl;
			return -2;
		}
		cerr << "[+] Wrote contigs to file (" << configuration.OutputFileName << ")." << endl;
		return 0;
	}
	return -1;
}
