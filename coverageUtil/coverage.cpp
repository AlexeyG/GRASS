/* 
 * File:   coverage.cpp
 * Author: alexeyg
 *
 * Created on December 3, 2011, 1:02 PM
 */

#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>
#include "Configuration.h"
#include "ReadCoverageReader.h"
#include "Reader.h"
#include "Sequence.h"

using namespace std;

typedef vector< int * > Depth;
typedef vector<FastASequence> Sequences;

Configuration config;
auto_ptr<ReadCoverage> coverage;
auto_ptr<Depth> depth;
auto_ptr<Sequences> contigs;


bool readContigs(const string &fileName, Sequences &contigs)
{
    FastAReader reader;
    bool result = reader.Open(fileName) && reader.Read(contigs) > 0;
    reader.Close();
    return result;
}

bool readCoverage(const string &fileName, ReadCoverage &coverage)
{
    ReadCoverageReader reader;
    bool result = reader.Open(fileName) && reader.Read(coverage);
    reader.Close();
    return result;
}

bool calculateDepth(const ReadCoverage &coverage, const Sequences &contigs, Depth &depth)
{
    int nContigs = coverage.GetContigCount();
    if (nContigs != (int)contigs.size())
        return false;
    depth.assign(nContigs, vector<int>());
    int avgReadLength = (int)coverage.AverageReadLength;
    for (int i = 0; i < nContigs; i++)
    {
        int contigLength = contigs[i].Nucleotides.length();
        depth[i].assign(contigLength, 0);
        for (vector<int>::const_iterator it = coverage.ReadLocations[i].begin(); it != coverage.ReadLocations[i].end(); it++)
        {
            for (int k = *it; k < *it + avgReadLength && k < contigLength; k++)
                depth[i][k]++;
        }
    }
    return true;
}

bool outputMIPSformat(const Sequences &contigs, const Depth &depth)
{
    int nContigs = contigs.size();
    for (int i = 0; i < nContigs; i++)
    {
        int contigLength = contigs[i].Nucleotides.length();
        long long coverageSum = 0;
        for (int j = 0; j < contigLength; j++)
            coverageSum += depth[i][j];
        double averageCoverage = (double)coverageSum / (double)contigLength;
        cout << i + 1 << "\t" << contigs[i].Comment << "\t" << contigLength << "\t" << averageCoverage << endl;
    }
    return true;
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(NULL));
    coverage = auto_ptr<ReadCoverage>(new ReadCoverage());
    depth = auto_ptr<Depth>(new Depth());
    contigs = auto_ptr<Sequences>(new Sequences());
    if (config.ProcessCommandLine(argc, argv))
    {
        if (!readContigs(config.ContigFileName, *contigs))
        {
            cerr << "[-] Unable to read contigs (" << config.ContigFileName << ")." << endl;
            return -2;
        }
        cerr << "[+] Read contigs (" << config.ContigFileName << ")." << endl;
        if (!readCoverage(config.CoverageFileName, *coverage))
        {
            cerr << "[-] Unable to read coverage (" << config.CoverageFileName << ")." << endl;
            return -3;
        }
        cerr << "[+] Read coverage (" << config.CoverageFileName << ")." << endl;
        if (!calculateDepth(*coverage, *contigs, *depth))
        {
            //clearDepthVector(depth);
            cerr << "[-] Unable to calculate coverage depth." << endl;
            return -4;
        }
        cerr << "[+] Calculated coverage depth." << endl;
        if (!outputMIPSformat(*contigs, *depth))
        {
            //clearDepthVector(depth);
            cerr << "[-] Unable to output coverage statistics." << endl;
            return -4;
        }
        //clearDepthVector(depth);
        return 0;
    }
    cerr << config.LastError;
    return -1;
}

