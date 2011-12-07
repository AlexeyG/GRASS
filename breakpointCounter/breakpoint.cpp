/* 
 * File:   breakpoint.cpp
 * Author: alexeyg
 *
 * Created on December 5, 2011, 7:42 PM
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <memory>
#include <boost/shared_array.hpp>
#include "Configuration.h"
#include "Sequence.h"
#include "Reader.h"
#include "Aligner.h"
#include "MummerCoord.h"
#include "MummerCoordReader.h"
#include "BreakpointCount.h"

using namespace std;

typedef vector<FastASequence> Sequences;
typedef vector<MummerCoord> Coords;

Configuration config;
auto_ptr<Sequences> scaffolds;
auto_ptr<Sequences> references;
auto_ptr<Coords> coords;
BreakpointCount breakpoints;


bool readContigs(const string &fileName, Sequences &contigs)
{
    FastAReader reader;
    bool result = reader.Open(fileName) && reader.Read(contigs) > 0;
    reader.Close();
    return result;
}

bool alignScaffolds(const string &referenceFileName, const string &scaffoldsFileName, const Sequences &references, const Sequences &scaffolds, Coords &coords)
{
    MummerAligner aligner(referenceFileName, scaffoldsFileName, config.MummerConfig);
    MummerCoordReader reader(references, scaffolds);
    if (!aligner.Align())
    {
        cerr << "[-] Unable to align scaffolds to reference (" << scaffoldsFileName << " -> " << referenceFileName << ")." << endl;
        return false;
    }
    aligner.RemoveOutput = false;
    if (!reader.Open(aligner.OutputFileName) || reader.Read(coords) == 0)
    {
        cerr << "[-] Unable to read MUMMER alignment (" << aligner.OutputFileName << ")." << endl;
        return false;
    }
    return true;
}

int filterAlignments(Coords &coords, double minBases)
{
    int count = 0;
    Coords newCoords;
    for (auto it = coords.begin(); it != coords.end(); it++)
        if (it->Identity * it->ReferenceAlignmentLength > minBases)
            newCoords.push_back(*it);
        else
            count++;
    coords.swap(newCoords);
    return count;
}

void output(const Coords &coords)
{
    for (Coords::const_iterator it = coords.begin(); it != coords.end(); it++)
        cout << it->QueryID << " [" << it->QueryPosition << ", " << it->QueryPosition + it->QueryAlignmentLength << "]" << endl;
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(NULL));
    scaffolds = auto_ptr<Sequences>(new Sequences());
    references = auto_ptr<Sequences>(new Sequences());
    coords = auto_ptr<Coords>(new Coords());
    if (config.ProcessCommandLine(argc, argv))
    {
        if (!readContigs(config.ScaffoldFileName, *scaffolds))
        {
            cerr << "[-] Unable to read scaffolds (" << config.ScaffoldFileName << ")." << endl;
            return -2;
        }
        cerr << "[+] Read scaffolds (" << config.ScaffoldFileName << ")." << endl;
        if (!readContigs(config.ReferenceFileName, *references))
        {
            cerr << "[-] Unable to read references (" << config.ReferenceFileName << ")." << endl;
            return -3;
        }
        cerr << "[+] Read references (" << config.ReferenceFileName << ")." << endl;
        if (!alignScaffolds(config.ReferenceFileName, config.ScaffoldFileName, *references, *scaffolds, *coords))
            return -4;
        BreakpointCount::Sort(*coords);
        cerr << "[+] Sorted MUMMER alignments." << endl;
        breakpoints.DistanceThreshold = config.DistanceThreshold;
        cerr << "[+] Aligned scaffolds to reference (" << config.ScaffoldFileName << " -> " << config.ReferenceFileName << ")." << endl;
        cerr << "[i] Filtered out " << filterAlignments(*coords, config.MinBases) << " alignments." << endl;
        cerr << "[i] Found a total of " << breakpoints.ProcessAlignments(*coords, *references, *scaffolds) << " breakpoints:" << endl;
        cerr << "    [i] Joins:\t" << breakpoints.Joins << endl;
        cerr << "    [i] Order:\t" << breakpoints.Order << endl;
        cerr << "    [i] Orientation:\t" << breakpoints.Orientation << endl;
        cerr << "    [i] Distance:\t" << breakpoints.Distance << endl;
        fprintf(stderr, "[i] Reference coverage:\t%.2lf %%\n", breakpoints.GetReferenceCoverage() * 100);
        fprintf(stderr, "[i] Scaffold coverage:\t%.2lf %%\n", breakpoints.GetScaffoldCoverage() * 100);
        return 0;
    }
    cerr << config.LastError;
    return -1;
}
