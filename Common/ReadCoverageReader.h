#ifndef _READCOVERAGEREADER_H
#define _READCOVERAGEREADER_H

#include <string>
#include <iostream>
#include "ReadCoverage.h"

using namespace std;

class ReadCoverageReader
{
public:
    ReadCoverageReader();
    virtual ~ReadCoverageReader();

public:
    bool Open(const string &fileName);
    bool Close();
    bool Read(ReadCoverage &store);

private:
    bool readHeader(int &nContigs, ReadCoverage &coverage);
    bool readContigs(int nContigs, ReadCoverage &coverage);
    bool readContig(int id, ReadCoverage &coverage);
        
protected:
    fstream in;
};
#endif
