/*
 * breakpointCounter : counts scaffold breakpoint using MUMmer local alignment.
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

void banner()
{
    cerr << "This program comes with ABSOLUTELY NO WARRANTY; see LICENSE for details." << endl;
    cerr << "This is free software, and you are welcome to redistribute it" << endl;
    cerr << "under certain conditions; see LICENSE for details." << endl;
    cerr << endl;
}

int main(int argc, char* argv[])
{
    banner();
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
        cout << "[+] Aligned scaffolds to reference (" << config.ScaffoldFileName << " -> " << config.ReferenceFileName << ")." << endl;
        cout << "[i] Filtered out " << filterAlignments(*coords, config.MinBases) << " alignments." << endl;
        cout << "[i] Found a total of " << breakpoints.ProcessAlignments(*coords, *references, *scaffolds) << " breakpoints:" << endl;
        cout << "    [i] Joins:\t" << breakpoints.Joins << endl;
        cout << "    [i] Order:\t" << breakpoints.Order << endl;
        cout << "    [i] Orientation:\t" << breakpoints.Orientation << endl;
        cout << "    [i] Distance:\t" << breakpoints.Distance << endl;
        fprintf(stdout, "[i] Reference coverage:\t%.2lf %%\n", breakpoints.GetReferenceCoverage() * 100);
        fprintf(stdout, "[i] Scaffold coverage:\t%.2lf %%\n", breakpoints.GetScaffoldCoverage() * 100);
        return 0;
    }
    cerr << config.LastError;
    return -1;
}
