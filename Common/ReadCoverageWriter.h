/* 
 * File:   ReadCoverageWriter.h
 * Author: alexeyg
 *
 * Created on November 29, 2011, 4:17 PM
 */

#ifndef _READCOVERAGEWRITER_H
#define	_READCOVERAGEWRITER_H

#include <string>
#include "ReadCoverage.h"

using namespace std;

class ReadCoverageWriter {
public:
    ReadCoverageWriter();
    virtual ~ReadCoverageWriter();

public:
    bool Open(const string &fileName, const string &mode = "wb");
    bool Close();
    bool Write(const ReadCoverage &coverage);

protected:
	FILE *out;
};

#endif	/* _READCOVERAGEWRITER_H */

