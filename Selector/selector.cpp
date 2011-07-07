/*
 * Data selector
 * Selects a subset of data (sequence, aligned reads, etc), which belong to a chosen segment of the genome.
 */

#include <iostream>
#include <cstdio>
#include <vector>
#include <set>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include "globals.h"
#include "configuration.h"
#include "helpers.h"
#include "reader.h"
#include "sequence.h"
#include "writer.h"
#include "XATag.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

#define VERSION "0.002"
#define AUTHOR "Alexey Gritsenko"
#define DATE "17.03.2011"

Configuration config;
vector<FastASequence> contigs;
set<Segment> segments;
vector<FastASequence> segmentSequence;

void printHelpMessage(void)
{
	cerr << "[i] Data selector version " << VERSION << " (" << DATE << ")" << endl;
	cerr << "[i] By " << AUTHOR << endl;
	cerr << "[i] Usage: dataSelector [arguments] <sequence.fasta>" << endl;
	cerr << "[i] -help                                               Print this message and exit." << endl;
	cerr << "[i] -chromosomes                                        Print chromosome numbers and exit." << endl;
	cerr << "[i] -segment <start> <length> [chromosome #]            Selects segment [start; start + length) from given chromosome. [1]" << endl;
	cerr << "                                                        Indices are 1-based." << endl;
	cerr << "[i] -select <length> <chromosome #> [count]             Randomly selects [count] regions of given length in the given chromosome." << endl;
	cerr << "[i] -selectall <length> [count]                         Randomly selects [count] regions of given length in every chromosome." << endl;
	cerr << "[i] -paired <input1.bam> <input2.bam> <output prefix>   Processes paired end BAM alignments and outputs two read files with given prefix." << endl;
	cerr << "[i] -output [sequence.fasta]                            Output filename for the selected sequences. [out.fasta]" << endl;
}

bool processCommandLine(int argc, char *argv[])
{
	if (argc == 1)
	{
		cerr << "[-] Not enough arguments. Consult -help." << endl;
		return false;
	}

	config.OutputFastaFileName = "out.fasta";

	int i = 1;
	while (i < argc)
	{
		if (!strcmp("-help", argv[i]) || !strcmp("-h", argv[i]))
		{
			printHelpMessage();
			return false;
		}
		else if (!strcmp("-chromosomes", argv[i]))
			config.PrintChromosomeInfo = true;
		else if (!strcmp("-segment", argv[i]))
		{
			if (argc - i - 1 < 2)
			{
				cerr << "[-] Parsing error in -segment: must have at least two arguments." << endl;
				return false;
			}
			i++;
			int start = atoi(argv[i]);
			if (start <= 0)
			{
				cerr << "[-] Parsing error in -segment: start must be a positive number." << endl;
				return false;
			}
			i++;
			int length = atoi(argv[i]); 
			if (length <= 0)
			{
				cerr << "[-] Parsing error in -segment: length must be a positive number." << endl;
				return false;
			}
			int chromosome = 1;
			if (i + 1 < argc && Helpers::IsNumber(argv[i + 1]))
			{
				i++;
				chromosome = atoi(argv[i]);
			}
			if (chromosome <= 0)
			{
				cerr << "[-] Parsing error in -segment: chromosome must be a positive number." << endl;
				return false;
			}
			start--;
			config.Segments.push_back(Segment(start, start + length - 1, chromosome));
		}
		else if (!strcmp("-select", argv[i]))
		{
			if (argc - i - 1 < 2)
			{
				cerr << "[-] Parsing error in -select: must have at least two argument." << endl;
				return false;
			}
			i++;
			int selectCount = 1;
			int selectChromosome = atoi(argv[i + 1]);
			int selectLength = atoi(argv[i]);
			if (selectLength <= 0)
			{
				cerr << "[-] Parsing error in -select: length must be a positive number." << endl;
				return false;
			}
			if (selectChromosome <= 0)
			{
				cerr << "[-] Parsing error in -select: chromosome # must be a positive number." << endl;
				return false;
			}
			i++;
			if (i + 1 < argc && Helpers::IsNumber(argv[i + 1]))
			{
				i++;
				selectCount = atoi(argv[i]);
			}
			if (selectCount <= 0)
			{
				cerr << "[-] Parsing error in -select: count must be a positive number." << endl;
				return false;
			}
			config.Select.push_back(ConfigSelect(selectLength, selectChromosome, selectCount));
		}
		else if (!strcmp("-selectall", argv[i]))
		{
			if (argc - i - 1 < 1)
			{
				cerr << "[-] Parsing error in -selectall: must have at least one argument." << endl;
				return false;
			}
			i++;
			int selectCount = 1;
			int selectChromosome = -1;
			int selectLength = atoi(argv[i]);
			if (selectLength <= 0)
			{
				cerr << "[-] Parsing error in -selectall: length must be a positive number." << endl;
				return false;
			}
			if (i + 1 < argc && Helpers::IsNumber(argv[i + 1]))
			{
				i++;
				selectCount = atoi(argv[i]);
			}
			if (selectCount <= 0)
			{
				cerr << "[-] Parsing error in -selectall: count must be a positive number." << endl;
				return false;
			}
			config.Select.push_back(ConfigSelect(selectLength, selectChromosome, selectCount));
		}
		else if (!strcmp("-output", argv[i]))
		{
			if (argc - i - 1 < 1)
			{
				cerr << "[-] Parsing error in -output: must have an argument." << endl;
				return false;
			}
			i++;
			config.OutputFastaFileName = argv[i];
		}
		else if (!strcmp("-paired", argv[i]))
		{
			if (argc - i - 1 < 3)
			{
				cerr << "[-] Parsing error in -paired: must have three arguments." << endl;
				return false;
			}
			i++;
			string inputBam1 = argv[i];
			i++;
			string inputBam2 = argv[i];
			i++;
			string outputPrefix = argv[i];
			config.PairedAlignment.push_back(PairedBam(inputBam1, inputBam2, outputPrefix));
		}
		else if (i == argc - 1)
			config.InputFastaFileName = argv[argc - 1];
		else
		{
			cerr << "[i] Unknown argument: " << argv[i] << endl;
			return false;
		}
		i++;
	}

	if (config.InputFastaFileName == "")
	{
		cerr << "[-] No sequence file given." << endl;
		return false;
	}

	return true;
}

bool readContigs(void)
{
	FastAReader reader;
	if (!reader.Open(config.InputFastaFileName))
		return false;

	if (reader.Read(contigs) == 0)
		return false;

	return true;
}

void printContigInformation()
{
	cerr << "[i] Contig information:" << endl;
	fprintf(stderr, "   [i] %2s %100s %15s\n", "#", "Name", "Length, bp");
	for (int i = 0; i < (int)contigs.size(); i++)
		fprintf(stderr, "   [i] %2i %100s %15i\n", i + 1, contigs[i].Comment.c_str(), (int)contigs[i].Nucleotides.length());
}

bool checkConfig(void)
{
	bool res = true;
	int nContigs = contigs.size();
	int nSegments = config.Segments.size();
	int nSelect = config.Select.size();

	cerr << "[i] Input data consistency check:" << endl;
	if (nSegments == 0 && nSelect == 0)
	{
		res = false;
		cerr << "   [-] Nothing to select." << endl;
	}

	for (int i = 0; i < nSegments; i++)
	{
		Segment seg = config.Segments[i];
		if (seg.Chromosome > nContigs)
		{
			res = false;
			cerr << "   [-] Segment " << i + 1 << ": belongs to a non-existent chromosome." << endl;
		}
		else
		{
			int len = contigs[seg.Chromosome - 1].Nucleotides.length();
			if (seg.Start > len || seg.Finish > len)
			{
				res = false;
				cerr << "   [-] Segment " << i + 1 << ": is out of chromosome boundaries." << endl;
			}
		}
	}

	for (int i = 0; i < nSelect; i++)
	{
		ConfigSelect sel = config.Select[i];
		if (sel.Chromosome > nContigs)
		{
			res = false;
			cerr << "   [-] Select " << i + 1 << ": belongs to a non-existent chromosome." << endl;
		}
	}

	vector<long long> selectLength(nContigs, 0);
	for (int i = 0; i < nSelect; i++)
	{
		ConfigSelect sel = config.Select[i];
		if (sel.Chromosome > nContigs)
			continue;
		if (sel.Chromosome > 0)
			selectLength[sel.Chromosome - 1] += sel.Length * (long long)sel.Count;
		else
			for (int j = 0; j < nContigs; j++)
				selectLength[j] += sel.Length;
	}
	for (int i = 0; i < nContigs; i++)
		if (selectLength[i] > (int)contigs[i].Nucleotides.length())
		{
			res = false;
			cerr << "   [-] Combined selection length for chromosome " << i + 1 << " is too large." << endl;
		}

	if (res)
		cerr << "   [+] Passed!" << endl;
	return res;
}

bool generateSegment(int length, int chromosome, Segment &seg)
{
	vector< pair<int, int> > usable;
	int last = 0;
	for (set<Segment>::iterator it = segments.begin(); it != segments.end(); it++)
	{
		if (it->Chromosome < chromosome)
			continue;
		if (it->Chromosome > chromosome)
			break;
		int start = it->Start;
		int finish = it->Finish;
		if (start - last - 1 >= length)
			usable.push_back(pair<int, int>(last, start - 1));
		last = finish + 1;
	}
	int chrLen = contigs[chromosome].Nucleotides.size();
	if (chrLen - last - 1 >= length)
		usable.push_back(pair<int, int>(last, chrLen));

	int count = usable.size();
	if (count == 0)
		return false;
	int id = rand() % count;
	
	seg.Start = rand() % (usable[id].second - usable[id].first + 1 - length) + usable[id].first;
	seg.Finish = seg.Start + length - 1;
	seg.Chromosome = chromosome;
	return true;
}

bool generateSegments(void)
{
	int nSegments = config.Segments.size();
	for (int i = 0; i < nSegments; i++)
	{
		Segment seg = config.Segments[i];
		segments.insert(Segment(seg.Start - 1, seg.Finish - 1, seg.Chromosome - 1));
	}

	Segment seg;
	int nSelect = config.Select.size();
	int nContigs = contigs.size();
	cerr << "[i] Generating segments:" << endl;
	for (int i = 0; i < nSelect; i++)
	{
		ConfigSelect sel = config.Select[i];
		if (sel.Chromosome > 0)
		{
			for (int i = 0; i < sel.Count; i++)
				if (!generateSegment(sel.Length, sel.Chromosome - 1, seg))
				{
					cerr << "   [-] Unable to add a segment of length " << sel.Length << " to chromosome " << sel.Chromosome << endl;
					return false;
				}
				else
				{
					segments.insert(seg);
					cerr << "   [+] Added segment [" << seg.Start << "; " << seg.Finish << "] to chromosome " << seg.Chromosome << endl;
				}
		}
		else
			for (int i = 0; i < nContigs; i++)
			{
				for (int j = 0; j < sel.Count; j++)
					if (!generateSegment(sel.Length, i, seg))
					{
						cerr << "   [-] Unable to add a segment of length " << sel.Length << " to chromosome " << sel.Chromosome << endl;
						return false;
					}
					else
					{
						segments.insert(seg);
						cerr << "   [+] Added segment [" << seg.Start << "; " << seg.Finish << "] to chromosome " << seg.Chromosome << endl;
					}
			}
	}
	return true;
}

void generateSegmentSequences()
{
	segmentSequence.resize(segments.size());
	char *buf = new char[MaxLine];
	int i = 0;
	for (set<Segment>::iterator it = segments.begin(); it != segments.end(); it++, i++)
	{
		Segment seg = *it;
		sprintf(buf, "%i-%i|%s", seg.Start, seg.Finish, contigs[seg.Chromosome].Comment.c_str());
		segmentSequence[i] = FastASequence(contigs[seg.Chromosome].Nucleotides.substr(seg.Start, seg.Finish - seg.Start + 1), buf);
	}
	delete [] buf;
}

bool outputSegmentSequences()
{
	FastAWriter writer;
	writer.Open(config.OutputFastaFileName);
	bool res = writer.Write(segmentSequence);
	writer.Close();

	return res;
}

void getReferenceIDs(BamReader &reader, vector<int> &ids)
{
	int nContigs = contigs.size();
	ids.resize(nContigs, -1);

	RefVector bam1vector = reader.GetReferenceData();
	for (RefVector::iterator it = bam1vector.begin(); it != bam1vector.end(); it++)
		for (int i = 0; i < nContigs; i++)
			if (it->RefName == contigs[i].Name())
			{
				ids[i] = reader.GetReferenceID(it->RefName);
				break;
			}
}

void convertToTags(const BamAlignment &alg, const BamReader &reader, vector<XATag> &tags)
{
	if (alg.IsMapped())
	{
		string xa;
		tags.push_back(XATag(alg));
		alg.GetTag("XA", xa);

		int i = 0, s = 0;
		int len = xa.length();
		while (i < len)
		{
			for (int t = 0; t < 3; t++)
			{
				while (i < len && xa[i] != ',') i++;
				i++;
			}
			while (i < len && xa[i] != ';') i++;
			if (i == len)
				break;
			i++;
			tags.push_back(XATag(xa.substr(s, s - i), reader));
			s = ++i;
		}
	}
}

vector<XATag> convertToTags(const BamAlignment &alg, const BamReader &reader)
{
	vector<XATag> tags;
	convertToTags(alg, reader, tags);
	return tags;
}

vector<XATag> convertToTags(const vector<BamAlignment> &alg, const BamReader &reader)
{
	vector<XATag> tags;
	for (vector<BamAlignment>::const_iterator it = alg.begin(); it != alg.end(); it++)
		convertToTags(*it, reader, tags);
	return tags;
}

bool alignmentOverlaps(const vector<XATag> tags, const vector<int> conv)
{
	if (tags.size() > MaxHits)
		return false;

	int n = tags.size();
	for (set<Segment>::iterator it = segments.begin(); it != segments.end(); it++)
		for (int i = 0; i < n; i++)
			if (conv[it->Chromosome] == tags[i].RefID && it->Start <= tags[i].Position + 1 && tags[i].Position + 1 <= it->Finish)
				return true;
	return false;
}

bool getNextAlignmentGroup(BamReader &bam, BamAlignment &buffer, vector<BamAlignment> &alg, bool core = false)
{
	int n = 0;
	alg.clear();
	if (!buffer.Name.empty())
	{
		alg.push_back(buffer), n++;
		buffer.Name.clear();
	}

	while (true)
	{
		 // EOF
		if (core)
		{
			if (!bam.GetNextAlignmentCore(buffer))
			{
				buffer.Name.clear();
				break;
			}
		}
		else
		{
			if (!bam.GetNextAlignment(buffer))
			{
				buffer.Name.clear();
				break;
			}
		}
		
		if (n == 0 || alg[n - 1].Name == buffer.Name)
		{
			alg.push_back(buffer), n++;
			buffer.Name.clear();
		}
		else
			break;
	}

	return alg.size() != 0;
}

bool processPairedAlignment(const PairedBam &bam)
{
	BamReader bam1, bam2;
	if (!bam1.Open(bam.InputBam1) || !bam2.Open(bam.InputBam2))
		return false;

	vector<int> bam1ref, bam2ref;
	getReferenceIDs(bam1, bam1ref);
	getReferenceIDs(bam2, bam2ref);
	
	vector<BamAlignment> alg1, alg2;
	BamAlignment buffer1, buffer2;
	buffer1.Name.clear();
	buffer2.Name.clear();
	bool error = false;
	while (!error)
	{
		bool read1 = getNextAlignmentGroup(bam1, buffer1, alg1, true);
		bool read2 = getNextAlignmentGroup(bam2, buffer2, alg2, true);
		if (!read1 && !read2)
			break;
		if (read1 ^ read2)
		{
			error = true;
			break;
		}
	}

	if (!error)
	{
		bam1.Rewind();	bam2.Rewind();
		buffer1.Name = buffer2.Name = "";
		FastQWriter w1, w2;
		bool success = w1.Open(bam.OutputPrefix + "_1.fastq") && w2.Open(bam.OutputPrefix + "_2.fastq");
		if (success)
		{
			while (true)
			{
				bool read1 = getNextAlignmentGroup(bam1, buffer1, alg1);
				bool read2 = getNextAlignmentGroup(bam2, buffer2, alg2);
				if (!read1 && !read2)
					break;
				string str1, str2;
				if (alignmentOverlaps(convertToTags(alg1, bam1), bam1ref) || alignmentOverlaps(convertToTags(alg2, bam2), bam2ref))
				{
					w1.Write(FastQSequence(alg1[0]));
					w2.Write(FastQSequence(alg2[0]));
				}
			}
		}
		w1.Close();
		w2.Close();
	}
	bam1.Close();
	bam2.Close();
	return !error;
}

bool processPairedAlignments()
{
	int n = config.PairedAlignment.size();
	for (int i = 0; i < n; i++)
		if (!processPairedAlignment(config.PairedAlignment[i]))
		{
			cerr << "   [-] Unable to process paired BAM with prefix " << config.PairedAlignment[i].OutputPrefix << endl;
			return false;
		}
		else
			cerr << "   [+] Processed paired BAM with prefix " << config.PairedAlignment[i].OutputPrefix << endl;
	return true;
}

int main(int argc, char *argv[])
{
	srand((unsigned int)time(NULL));
	if (processCommandLine(argc, argv))
	{
		if (!readContigs())
		{
			cerr << "[-] Unable to open sequence file " << config.InputFastaFileName << endl;
			return -1;
		}
		if (config.PrintChromosomeInfo)
		{
			printContigInformation();
			return 0;
		}
		if (!checkConfig())
			return -2;
		if (!generateSegments())
			return -3;
		generateSegmentSequences();
		cerr << "[i] Outputting segment sequences to file " << config.OutputFastaFileName << endl;
		if (!outputSegmentSequences())
		{
			cerr << "   [-] Unable to write segment sequences." << endl;
			return -4;
		}
		cerr << "[i] Processing paired alignments." << endl;
		if (!processPairedAlignments())
			return -5;
		return 0;
	}
	return -1;
}
