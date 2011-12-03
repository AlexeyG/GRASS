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

Configuration config;
ReadCoverage coverage;
vector<FastASequence> contigs;
vector<int *> depth;


bool readContigs(const string &fileName, vector<FastASequence> &contigs)
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

bool calculateDepth(const ReadCoverage &coverage, const vector<FastASequence> &contigs, vector<int*> &depth)
{
    int nContigs = coverage.GetContigCount();
    if (nContigs != (int)contigs.size())
        return false;
    depth.resize(nContigs, NULL);
    int avgReadLength = (int)coverage.AverageReadLength;
    for (int i = 0; i < nContigs; i++)
    {
        int contigLength = contigs[i].Nucleotides.length();
        depth[i] = new int[contigLength];
        int readCount = (int)coverage.ReadLocations[i].size();
        for (int j = 0; j < readCount; j++)
        {
            int pos = coverage.ReadLocations[i][j];
            for (int k = pos; k < pos + avgReadLength && k < contigLength; k++)
                depth[i][k]++;
        }
    }
    return true;
}

bool outputMIPSformat(const vector<FastASequence> &contigs, const vector<int*> &depth)
{
    int nContigs = contigs.size();
    for (int i = 0; i < nContigs; i++)
    {
        int contigLength = contigs[i].Nucleotides.length();
        long long coverageSum = 0;
        for (int j = 0; j < contigLength; j++)
            coverageSum += depth[i][j];
        double averageCoverage = (double)coverageSum / (double)contigLength;
    }
    return true;
}

void clearDepthVector(vector<int *> &depth)
{
    int n = depth.size();
    for (int i = 0; i < n; i++)
        delete depth[i];
    depth.clear();
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(NULL));
    if (config.ProcessCommandLine(argc, argv))
    {
        if (!readContigs(config.ContigFileName, contigs))
        {
            cerr << "[-] Unable to read contigs (" << config.ContigFileName << ")." << endl;
            return -2;
        }
        cerr << "[+] Read contigs (" << config.ContigFileName << ")." << endl;
        if (!readCoverage(config.CoverageFileName, coverage))
        {
            cerr << "[-] Unable to read coverage (" << config.CoverageFileName << ")." << endl;
            return -3;
        }
        cerr << "[+] Read coverage (" << config.CoverageFileName << ")." << endl;
        if (!calculateDepth(coverage, contigs, depth))
        {
            cerr << "[-] Unable to calculate coverage depth." << endl;
            clearDepthVector(depth);
            return -4;
        }
        cerr << "[+] Calculated coverage depth." << endl;
        if (!outputMIPSformat(contigs, depth))
        {
            cerr << "[-] Unable to output coverage statistics." << endl;
            clearDepthVector(depth);
            return -4;
        }
        
        clearDepthVector(depth);
        return 0;
    }
    cerr << config.LastError;
    return -1;
}

