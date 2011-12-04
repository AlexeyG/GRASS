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
#include <boost/shared_array.hpp>
#include <cstring>
#include "Configuration.h"
#include "ReadCoverageReader.h"
#include "Reader.h"
#include "Sequence.h"

using namespace std;

typedef vector< boost::shared_array<int> > Depth;
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
    //depth.assign(nContigs, vector<int>());
    int avgReadLength = (int)coverage.AverageReadLength;
    for (int i = 0; i < nContigs; i++)
    {
        int contigLength = contigs[i].Nucleotides.length();
        depth[i] = boost::shared_array<int>(new int[contigLength]);
        memset(depth[i].get(), 0, contigLength * sizeof(int));
        for (vector<int>::const_iterator it = coverage.ReadLocations[i].begin(); it != coverage.ReadLocations[i].end(); it++)
        {
            for (int k = max(*it, 0); k < *it + avgReadLength && k < contigLength; k++)
                depth[i][k]++;
        }
    }
    return true;
}

bool outputDepth(const string &fileName, const Sequences &contigs, const Depth &depth)
{
    FILE *out = fopen(fileName.c_str(), "w");
    if (out == NULL)
        return false;
    
    int nContigs = contigs.size();
    for (int i = 0; i < nContigs; i++)
    {
        int contigLength = contigs[i].Nucleotides.length();
        for (int j = 0; j < contigLength; j++)
            fprintf(out, "%i\n", depth[i][j]);
    }
    fclose(out);
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
    contigs = auto_ptr<Sequences>(new Sequences());
    if (config.ProcessCommandLine(argc, argv))
    {
        if (!readContigs(config.ContigFileName, *contigs))
        {
            cerr << "[-] Unable to read contigs (" << config.ContigFileName << ")." << endl;
            return -2;
        }
        depth = auto_ptr<Depth>(new Depth(contigs->size(), boost::shared_array<int>())); // otherwise we get a segfault. Resizing vector in place in a pain (stack overrun?)
        cerr << "[+] Read contigs (" << config.ContigFileName << ")." << endl;
        if (!readCoverage(config.CoverageFileName, *coverage))
        {
            cerr << "[-] Unable to read coverage (" << config.CoverageFileName << ")." << endl;
            return -3;
        }
        cerr << "[+] Read coverage (" << config.CoverageFileName << ")." << endl;
        if (!calculateDepth(*coverage, *contigs, *depth))
        {
            cerr << "[-] Unable to calculate coverage depth." << endl;
            return -4;
        }
        cerr << "[+] Calculated coverage depth." << endl;
        if (!config.DepthFileName.empty())
        {
            if (!outputDepth(config.DepthFileName, *contigs, *depth))
            {
                cerr << "[-] Unable to output coverage depth (" << config.DepthFileName << ")." << endl;
                return -5;
            }
            cerr << "[+] Output coverage depth (" << config.DepthFileName << ")." << endl;
        }
        if (!outputMIPSformat(*contigs, *depth))
        {
            cerr << "[-] Unable to output coverage statistics." << endl;
            return -4;
        }
        return 0;
    }
    cerr << config.LastError;
    return -1;
}

