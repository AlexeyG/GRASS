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

#include <iostream>
#include <cstdio>
#include <vector>
#include <set>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include "Globals.h"
#include "Defines.h"
#include "Configuration.h"
#include "Reader.h"
#include "Sequence.h"
#include "Writer.h"
#include "XATag.h"
#include "Helpers.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

Configuration config;
vector<FastASequence> contigs;
set<Segment> segments;
vector<FastASequence> segmentSequence;

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

bool generatePairedReads(const PairedSimulation &simulation, int &counter)
{
	FastQWriter w1, w2;
	bool success = w1.Open(simulation.OutputPrefix + "_1.fastq") && w2.Open(simulation.OutputPrefix + "_2.fastq");
	if (success)
	{
		for (int i = 0; i < (int)segmentSequence.size(); i++)
		{
			int length = segmentSequence[i].Nucleotides.length();
			int readCount = (length * simulation.Depth) / (2 * simulation.ReadLengthMean);
			//printf("Mean %.2lf, Std %.2lf\n", simulation.ReadLengthMean, simulation.ReadLengthStd);
			while (readCount-- > 0)
			{
				int len1, len2, insert;
				len1 = len2 = insert = 0;
				while (len1 <= 0)
					len1 = Helpers::RandomNormal(simulation.ReadLengthMean, simulation.ReadLengthStd);
				while (len2 <= 0)
					len2 = Helpers::RandomNormal(simulation.ReadLengthMean, simulation.ReadLengthStd);
				while (insert <= 0)
					insert = Helpers::RandomNormal(simulation.InsertSizeMean, simulation.InsertSizeStd);
				if (insert < len1 + len2)
					insert = len1 + len2;
				int pos1 = rand() % length;
				///printf("%i : %i-%i %i %i-%i\n", length, pos1, len1, insert, pos1 + len1 + insert, len2);
				//printf("Length: %i --- %i\n", len1, len2);
				if (pos1 + insert < length)
				{
					FastQSequence l,r;
					string seq1 = segmentSequence[i].Nucleotides.substr(pos1, len1);
					string seq2 = segmentSequence[i].Nucleotides.substr(pos1 + insert - len2, len2);
					string qual1(len1, '~'), qual2(len2, '~');
					if (simulation.IsIllumina)
					{
						l = FastQSequence(seq1, Helpers::ItoStr(counter) + ".1|" + Helpers::ItoStr(pos1), qual1), r = FastQSequence(seq2, Helpers::ItoStr(counter) + ".2|" + Helpers::ItoStr(pos1 + insert - len2), qual2);
						r.ReverseCompelement();
					}
					else
						l = FastQSequence(seq2, Helpers::ItoStr(counter) + ".1|" + Helpers::ItoStr(pos1), qual2), r = FastQSequence(seq1, Helpers::ItoStr(counter) + ".2|" + Helpers::ItoStr(pos1 + insert - len2), qual1);
					
					/* Read flipping */
					if (rand() < RAND_MAX / 2)
					{
						FastQSequence t(l);
						l = r;
						r = t;

						if (!simulation.IsIllumina)
						{
							l.ReverseCompelement();
							r.ReverseCompelement();
						}
					}
					w1.Write(l), w2.Write(r);
					counter++;
				}
			}
		}
	}
	w1.Close();
	w2.Close();
	return success;
}

bool generatePairedReads()
{
	int n = config.PairedReadSimulation.size();
	int counter = 0;
	for (int i = 0; i < n; i++)
		if (!generatePairedReads(config.PairedReadSimulation[i], counter))
		{
			cerr << "   [-] Unable to generate paired reads with prefix " << config.PairedReadSimulation[i].OutputPrefix << endl;
			return false;
		}
		else
			cerr << "   [+] Generated paired reads with prefix " << config.PairedReadSimulation[i].OutputPrefix << endl;
	return true;
}

int main(int argc, char *argv[])
{
	srand((unsigned int)time(NULL));
	if (config.ProcessCommandLine(argc, argv))
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
		cerr << "[i] Generating simulated paired reads." << endl;
		if (!generatePairedReads())
			return -6;
		return 0;
	}
	cerr << config.LastError;
	return -1;
}
