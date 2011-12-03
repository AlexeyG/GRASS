/* 
 * File:   coverage.cpp
 * Author: alexeyg
 *
 * Created on December 3, 2011, 1:02 PM
 */

#include <cstdlib>
#include <iostream>
#include <vector>
#include "Configuration.h"
#include "ReadCoverageReader.h"
#include "Reader.h"
#include "Sequence.h"

using namespace std;

Configuration config;
ReadCoverage coverage;
vector<FastASequence> contigs;
vector< vector<int> > depth;


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

bool calculateDepth(const ReadCoverage &coverage, const vector<FastASequence> &contigs, vector< vector<int> > &depth)
{
    int nContigs = coverage.GetContigCount();
    if (nContigs != (int)contigs.size())
        return false;
    depth.assign(nContigs, vector<int>());
    cout << nContigs << " " << contigs.size() << " " << coverage.ReadLocations.size() << endl;
    int avgReadLength = (int)coverage.AverageReadLength;
    for (int i = 0; i < nContigs; i++)
    {
        cout << "A" << endl;
        int contigLength = contigs[i].Nucleotides.length();
        cout << "B" << endl;
        depth[i].assign(contigLength, 0);
        cout << "C" << endl;
        int readCount = (int)coverage.ReadLocations[i].size();
        cout << "D" << endl;
        cout << "Contig " << i << " of length " << contigLength << " with count " << readCount << endl;
        for (int j = 0; j < readCount; j++)
        {
            int pos = coverage.ReadLocations[i][j];
            //cout << "  Read " << j << " with position " << pos << " and length " << avgReadLength << endl;
            for (int k = pos; k < pos + avgReadLength && k < contigLength; k++)
                depth[i][k]++;
        }
    }
    return true;
}

bool outputMIPSformat(const vector<FastASequence> &contigs, const vector< vector<int> > &depth)
{
    int nContigs = contigs.size();
    for (int i = 0; i < nContigs; i++)
    {
        int contigLength = contigs[i].Nucleotides.length();
        long long coverageSum = 0;
        for (int j = 0; j < contigLength; j++)
            coverageSum += depth[i][j];
        double averageCoverage = (double)coverageSum / (double)contigLength;
        cout << i + i << "\t" << contigs[i].Comment << "\t" << contigLength << "\t" << averageCoverage << endl;
    }
    return true;
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
            return -4;
        }
        cerr << "[+] Calculated coverage depth." << endl;
        if (!outputMIPSformat(contigs, depth))
        {
            cerr << "[-] Unable to output coverage statistics." << endl;
            return -4;
        }
        
        return 0;
    }
    cerr << config.LastError;
    return -1;
}

