#include "MummerCoordReader.h"
#include "Globals.h"

#include <iostream>

using namespace std;

MummerCoordReader::MummerCoordReader()
{
    fin = NULL;
    line = new char[MaxLine];
    numCoords = -1;
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
    cout << "Want to open: " << fileName << endl;
    fin = fopen(fileName.c_str(), mode.c_str());
    cout << "Got: " << (int)fin << endl;
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
    cout << "Read: " << line << endl;
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
    fseek(fin, 0, SEEK_SET);
    
    if (numCoords >= 0)
        return numCoords;

    long long num = 0;
    while (fgets(line, MaxLine, fin) != NULL)
        ++num;

    fseek(fin, index, SEEK_SET);
    return numCoords = num;
}
