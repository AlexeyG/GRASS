/*
 * dataFilter : a support tool used in debuging to filter paired reads aligning
 * to gaps between contigs.
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
#include "ContigInfo.h"
#include "Reader.h"
#include "Writer.h"
#include "AlignmentReader.h"
#include <iostream>
#include <algorithm>

using namespace std;

Configuration config;
vector<ContigInfo> contigs;

bool readContigs(const string &fileName, vector<ContigInfo> &contigs)
{
	bool result = true;
	FastAReader reader;
	FastASequence seq;
	if (!reader.Open(fileName))
		result = false;
	while (reader.Read(seq))
		contigs.push_back(seq);
	reader.Close();
	sort(contigs.begin(), contigs.end());
	return result;
}

bool containedInContig(const BamAlignment &alg, const vector<XATag> &tags, const vector<ContigInfo> &contigs)
{
	vector<int> pos;
	pos.push_back(alg.Position);
	for (int i = 0; i < (int)tags.size(); i++)
		pos.push_back(tags[i].Position); 
	int n = contigs.size();
	int readLen = alg.Length;
	for (int j = 0; j < (int)pos.size(); j++)
	{
		bool contained = false;
		for (int i = 0; i < n; i++)
		{
			if (contigs[i].Position < 0)
				continue;
			if (pos[j] >= contigs[i].Position && pos[j] + readLen <= contigs[i].Position + contigs[i].Length)
				contained = true;
		}
		if (!contained)
			return false;
	}
	return true;
}

bool processPairedReads(const PairedInput &input, const vector<ContigInfo> &contigs)
{
	int readsInGaps = 0, pairsInGaps = 0;
	bool success = true;
	AlignmentReader left, right;
	FastQWriter leftOut, rightOut;
	if (!left.Open(input.LeftFileName) || !right.Open(input.RightFileName))
		success = false;
	if (success && (!leftOut.Open(input.OutputPrefix + "_1.fastq") || !rightOut.Open(input.OutputPrefix + "_2.fastq")))
		success = false;
	if (success)
	{
		vector<XATag> leftTags, rightTags;
		BamAlignment leftAlignment, rightAlignment;
		while (left.GetNextAlignmentGroup(leftAlignment, leftTags) && right.GetNextAlignmentGroup(rightAlignment, rightTags))
		{
			bool leftContained = containedInContig(leftAlignment, leftTags, contigs);
			bool rightContained = containedInContig(rightAlignment, rightTags, contigs);
			if (leftContained && rightContained)
			{
				FastQSequence leftSeq(leftAlignment), rightSeq(rightAlignment);
				leftOut.Write(leftSeq), rightOut.Write(rightSeq);
			}
			if (!leftContained)
				readsInGaps++;
			if (!rightContained)
				readsInGaps++;
			if (!leftContained || !rightContained)
				pairsInGaps++;
		}
	}
	left.Close();
	right.Close();
	leftOut.Close();
	rightOut.Close();
	cout << input.OutputPrefix << ": " << readsInGaps << " " << pairsInGaps << endl;
	return success;
}

bool processPairedReads(const vector<PairedInput> &input, const vector<ContigInfo> &contigs)
{
	int n = input.size();
	for (int i = 0; i < n; i++)
	{
		cerr << "    [i] Processing paired input: " << input[i].OutputPrefix << endl;
		if (!processPairedReads(input[i], contigs))
			return false;
	}
	return true;
}

int main(int argc, char *argv[])
{
	if (config.ProcessCommandLine(argc, argv))
	{
		if (!readContigs(config.InputFileName, contigs))
		{
			cerr << "[-] Unable to read contigs: " << config.InputFileName << endl;
			return -2;
		}
		cerr << "[+] Successfully read contigs." << endl;
		cerr << "[i] Processing paired reads:" << endl;
		if (!processPairedReads(config.PairedFilter, contigs))
		{
			cerr << "[-] Unable to process paired reads." << endl;
			return -3;
		}
		cerr << "[+] Successfully processed paired reads." << endl;
		return 0;
	}
	cerr << config.LastError;
	return -1;
}
