/*
 * kmer : creates k-mer count statistics for a given k.
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
#include "Reader.h"
#include "Configuration.h"
#include "Location.h"
#include <map>

using namespace std;

Configuration config;
vector<FastASequence> contigs;
map<string, int> kmers;
map<string, vector<Location> > locations;

bool readContigs(const string &filename, vector<FastASequence> &contigs)
{
	FastAReader reader;
	bool success = reader.Open(filename);
	int read = -1;
	if (success)
		read = reader.Read(contigs);
	reader.Close();
	return success && read > 0;
}

void countKmers(int k, int id, const FastASequence &seq, map<string, int> &kmers, map<string, vector<Location> > &locations)
{
	int len = seq.Nucleotides.length();
	for (int i = 0; i < len - k; i++)
	{
		string kmer = seq.Nucleotides.substr(i, k);
		if (kmers.find(kmer) == kmers.end())
		{
			kmers[kmer] = 1;
			locations[kmer] = vector<Location>();
		}
		else
			kmers[kmer]++;
		locations[kmer].push_back(Location(id, i));
	}
}

void countKmers(int k, const vector<FastASequence> &contigs, map<string, int> &kmers, map<string, vector<Location> > &locations)
{
	int n = contigs.size();
	for (int i = 0; i < n; i++)
	{
		cerr << "[i] Adding contig " << i + 1 << "." << endl;
		countKmers(k, i, contigs[i], kmers, locations);
	}
}

void printKmers(const map<string, int> &kmers, bool verbose)
{
	int num = 0;
	int count = 0;
	for (map<string, int>::const_iterator it = kmers.begin(); it != kmers.end(); it++)
		if (it->second > 1)
		{
			if (verbose)
			{
				printf("   [i] Duplicate: %s %i\n", it->first.c_str(), it->second);
				printf("   [i] ");
				for (vector<Location>::const_iterator loc = locations[it->first].begin(); loc != locations[it->first].end(); loc++)
					printf("%s ", loc->ToString().c_str());
				printf("\n");
			}
			count += it->second;
			num++;
		}
	printf("   [i] Duplicates: %i\n", num);
	printf("   [i] Copies: %i\n", count);
}

int main(int argc, char *argv[])
{
	if (config.ProcessCommandLine(argc, argv))
	{
		if (!readContigs(config.InputFileName, contigs))
		{
			cerr << "[-] Unable to read contigs: " << config.InputFileName << endl;
			return -1;
		}
		countKmers(config.Kmer, contigs, kmers, locations);
		cerr << "[+] Outputing stats:" << endl;
		printKmers(kmers, config.Verbose);
		return 0;
	}
	cerr << config.LastError;
	return -1;
}
