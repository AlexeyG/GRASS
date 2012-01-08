#ifndef _MUMMERTILINGREADER_H
#define	_MUMMERTILINGREADER_H

#include <vector>
#include <cstdio>
#include <string>
#include <map>
#include "MummerTiling.h"
#include "Sequence.h"

using namespace std;

class MummerTilingReader
{
public:
    MummerTilingReader(const vector<FastASequence> &references, const vector<FastASequence> &scaffolds);
    ~MummerTilingReader();
    
public:
    bool Open(const string &fileName, const string &mode = "rb");
    bool Close();
    bool IsOpen() const;
    bool Read(MummerTiling &tiling);
    long long Read(vector<MummerTiling> &tilings);
    long long NumTilings();
    
private:
    FILE *fin;
    char *line;
    long long numCoords;
    int referenceID;
    
private:
    map<string, int> referenceIds;
    map<string, int> scaffoldIds;
    
private:
    void createMap(map<string, int> &store, const vector<FastASequence> &seq);
};

#endif

