/*
 * Common : a collection of classes (re)used throughout the scaffolder implementation.
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


#include "MummerTilingReader.h"
#include "Globals.h"
#include <cstring>
#include <sstream>

using namespace std;

MummerTilingReader::MummerTilingReader(const vector<FastASequence> &references, const vector<FastASequence> &scaffolds)
{
    fin = NULL;
    line = new char[MaxLine];
    numCoords = -1;
    referenceID = -1;
    createMap(referenceIds, references);
    createMap(scaffoldIds, scaffolds);
}

MummerTilingReader::~MummerTilingReader()
{
    Close();
    delete [] line;
}

bool MummerTilingReader::Open(const string &fileName, const string &mode)
{
    if (fin != NULL)
        return false;
    fin = fopen(fileName.c_str(), mode.c_str());
    if (fin == NULL)
        return false;
    return true;
}

bool MummerTilingReader::Close()
{
    if (fin != NULL)
    {
        fclose(fin);
        fin = NULL;
        numCoords = -1;
        referenceID = -1;
        return true;
    }
    return false;
}

bool MummerTilingReader::IsOpen() const
{
    return fin != NULL;
}

bool MummerTilingReader::Read(MummerTiling &tiling)
{
    do
    {
        if (fgets(line, MaxLine, fin) == NULL)
            return false;

        int len = strlen(line);
        if (len <= 1)
            return false;
        line[--len] = '\0';

        if (line[0] == '>')
        {
            stringstream tmpStream(line + 1);
            string newReference;
            tmpStream >> newReference;
            auto it = referenceIds.find(newReference);
            referenceID = (it == referenceIds.end() ? -1 : it->second);
        }
    } while (line != NULL && line[0] == '>');
    
    // Columns are, start in ref, end in ref
    // distance to next contig, length of this contig, alignment
    // coverage, identity, orientation, and ID respectively.
    
    int t, endPosition;
    stringstream ss(line);
    tiling.ReferenceID = referenceID;
    ss >> tiling.ReferencePosition >> endPosition >> t >> tiling.QueryLength >> tiling.Coverage >> tiling.Identity;
    tiling.ReferenceLength = endPosition - tiling.ReferencePosition + 1;
    tiling.Coverage /= 100;
    tiling.Identity /= 100;
    
    string tStr, orientationString, queryName;
    getline(ss, tStr, '\t');
    getline(ss, orientationString, '\t');
    getline(ss, queryName, '\t');
    
    tiling.IsReverse = orientationString == "-";
    //cout << "Orientation: <" << orientationString << "> length = " << tiling.ReferenceLength << endl;
    auto it = scaffoldIds.find(queryName); tiling.QueryID = (it == scaffoldIds.end() ? -1 : it->second);
    
    // Turning to 0-based
    tiling.ReferencePosition--;
    
    return true;
}

long long MummerTilingReader::Read(vector<MummerTiling> &tilings)
{
    tilings.resize(NumTilings());
    for (unsigned i = 0; i < tilings.size(); ++i)
        Read(tilings[i]);

    return tilings.size();
}

long long MummerTilingReader::NumTilings()
{
    if (!IsOpen())
        return 0;
    long long index = ftell(fin);
    
    if (numCoords >= 0)
        return numCoords;

    fseek(fin, 0, SEEK_SET);
    long long num = 0;
    while (fgets(line, MaxLine, fin) != NULL)
        if (line[0] != '>')
            ++num;
    
    fseek(fin, index, SEEK_SET);
    return numCoords = num;
}

void MummerTilingReader::createMap(map<string, int> &store, const vector<FastASequence> &seq)
{
    int nSeq = seq.size();
    for (int i = 0; i < nSeq; i++)
        store[seq[i].Name()] = i;
}