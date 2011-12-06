#include "MummerCoordReader.h"
#include "Globals.h"
#include <sstream>

#include <iostream>

using namespace std;

MummerCoordReader::MummerCoordReader(const vector<FastASequence> &references, const vector<FastASequence> &scaffolds)
{
    fin = NULL;
    line = new char[MaxLine];
    numCoords = -1;
    createMap(referenceIds, references);
    createMap(scaffoldsIds, scaffolds);
}

MummerCoordReader::~MummerCoordReader()
{
    Close();
    delete [] line;
}

bool MummerCoordReader::Open(const string &fileName, const string &mode)
{
    if (fin != NULL)
        return false;
    fin = fopen(fileName.c_str(), mode.c_str());
    if (fin == NULL)
        return false;
    return true;
}

bool MummerCoordReader::Close()
{
    if (fin != NULL)
    {
        fclose(fin);
        fin = NULL;
        numCoords = -1;
        return true;
    }
    return false;
}

bool MummerCoordReader::IsOpen() const
{
    return fin != NULL;
}

bool MummerCoordReader::Read(MummerCoord &coord)
{
    if (fgets(line, MaxLine, fin) == NULL)
        return false;
    
    int t;
    stringstream ss(line);
    ss >> coord.ReferencePosition >> t >> coord.QueryPosition >> t >> coord.ReferenceAlignmentLength >> coord.QueryAlignmentLength >> coord.Identity;
    coord.Identity /= 100;
    ss >> t;
    coord.IsReferenceReverse = t > 0;
    ss >> t;
    coord.IsQueryReverse = t > 0;
    
    cout << "LEFT WITH:" << ss.str() << endl;
    
    return true;
}

long long MummerCoordReader::Read(vector<MummerCoord> &coords)
{
    coords.resize(NumCoords());
    for (unsigned i = 0; i < coords.size(); ++i)
        Read(coords[i]);

    return coords.size();
}

long long MummerCoordReader::NumCoords()
{
    if (!IsOpen())
        return 0;
    long long index = ftell(fin);
    
    if (numCoords >= 0)
        return numCoords;

    fseek(fin, 0, SEEK_SET);
    long long num = 0;
    while (fgets(line, MaxLine, fin) != NULL)
        ++num;
    
    fseek(fin, index, SEEK_SET);
    return numCoords = num;
}

void createMap(map<string, int> &store, const vector<FastASequence> &seq)
{
    int nSeq = seq.size();
    for (int i = 0; i < nSeq; i++)
        store[seq[i].Name()] = i;
}
