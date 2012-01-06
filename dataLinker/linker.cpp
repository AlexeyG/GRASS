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

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include "Configuration.h"
#include "DataStore.h"
#include "PairedReadConverter.h"
#include "SequenceConverter.h"
#include "DataStoreWriter.h"
#include "ReadCoverage.h"
#include "ReadCoverageWriter.h"

using namespace std;

Configuration config;
DataStore store;
ReadCoverage coverage;

bool processPairs(const Configuration &config, DataStore &store, const vector<PairedInput> &paired, ReadCoverage &coverage)
{
	PairedReadConverter converter(store);
	int n = (int)paired.size();
	for (int i = 0; i < n; i++)
	{
		PairedInput p = paired[i];
		cerr << "   [i] Processing " << (p.IsIllumina ? "Illumina" : "454") << " paired reads (" << p.LeftFileName << ", " << p.RightFileName << ") of weight " << p.Weight << " with insert size " << p.Mean << " +/- " << p.Std << endl; 
		switch (converter.Process(config, paired[i]))
		{
		case PairedReadConverter::Success:
			cerr << "      [+] Successfully processed paired reads." << endl;
			break;
		case PairedReadConverter::FailedLeftAlignment:
			cerr << "      [-] Unable to align left read mates (" << (p.IsIllumina ? "Illumina" : "454") << ")." << endl;
			return false;
		case PairedReadConverter::FailedRightAlignment:
			cerr << "      [-] Unable to align second read mates (" << (p.IsIllumina ? "Illumina" : "454") << ")." << endl;
			return false;
		case PairedReadConverter::FailedLeftConversion:
			cerr << "      [-] Unable to convert left alignment SAM to BAM (SAM Tools)." << endl;
			return false;
		case PairedReadConverter::FailedRightConversion:
			cerr << "      [-] Unable to convert right alignment SAM to BAM (SAM Tools)." << endl;
			return false;
		case PairedReadConverter::FailedLinkCreation:
			cerr << "      [-] Unable to create contig links from alignment." << endl;
			return false;
                case PairedReadConverter::InconsistentReferenceSets:
                        cerr << "      [-] Inconsistent reference sets in alignments." << endl;
                        return false;
		}
	}
        coverage = converter.ContigReadCoverage;
	return true;
}

bool processSequences(const Configuration &config, DataStore &store, const vector<SequenceInput> &sequences)
{
    SequenceConverter converter(store);
    int n = (int)sequences.size();
    for (int i = 0; i < n; i++)
    {
        SequenceInput s = sequences[i];
        cerr << "   [i] Processing sequences (" << s.FileName << ") of weight " << s.Weight << " and deviation " << s.Std << endl;
        switch (converter.Process(config, s))
        {
            case SequenceConverter::Success:
                cerr << "      [+] Successfully processed sequences." << endl;
                break;
            case SequenceConverter::FailedAlignment:
                cerr << "      [-] Unable to align contigs." << endl;
                return false;
            case SequenceConverter::FailedReadAlignment:
                cerr << "      [-] Unable to read contig alignment." << endl;
                return false;
            case SequenceConverter::FailedReadSequences:
                cerr << "      [-] Unable to read sequences." << endl;
                return false;
        }
    }
    return true;
}

bool writeStore(const DataStore &store, const string &fileName)
{
	DataStoreWriter writer;
	bool result = writer.Open(fileName) && writer.Write(store);
	writer.Close();
	return result;
}

bool writeCoverage(const ReadCoverage &coverage, const string &fileName)
{
    ReadCoverageWriter writer;
    bool result = writer.Open(fileName) && writer.Write(coverage);
    writer.Close();
    return result;
}

int main(int argc, char *argv[])
{
	srand((unsigned int)time(NULL));
	if (config.ProcessCommandLine(argc, argv))
	{
		if (!store.ReadContigs(config.InputFileName))
		{
                    cerr << "[-] Unable to read contigs from file (" << config.InputFileName << ")." << endl;
                    return -2;
		}
		cerr << "[+] Read input contigs from file (" << config.InputFileName << ")." << endl;
		cerr << "[i] Processing paired reads." << endl;
		if (!processPairs(config, store, config.PairedReadInputs, coverage))
			return -3;
                if (!processSequences(config, store, config.SequenceInputs))
			return -4;
		if (!writeStore(store, config.OutputFileName))
		{
                    cerr << "[-] Unable to output generated links into file (" << config.OutputFileName << ")." << endl;
                    return -5;
		}
                cerr << "[+] Output generated link into file (" << config.OutputFileName << ")." << endl;
                if (!config.ReadCoverageFileName.empty())
                {
                    if (!writeCoverage(coverage, config.ReadCoverageFileName))
                    {
                        cerr << "[-] Unable to output read coverage statistics into file (" << config.ReadCoverageFileName << ")." << endl;
                        return -6;
                    }
                    cerr << "[+] Output read coverage into file (" << config.ReadCoverageFileName << ")." << endl;
                }
		return 0;
	}
	cerr << config.LastError;
	return -1;
}
