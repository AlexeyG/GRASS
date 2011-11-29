#ifndef _ALIGNMENTREADER_H
#define _ALIGNMENTREADER_H

#include <iostream>
#include "api/BamReader.h"
#include "XATag.h"

using namespace std;
using namespace BamTools;

class AlignmentReader
{
public:
	AlignmentReader();

public:
	bool Open(const string &fileName);
	bool Close();
	bool IsOpen() const;
	bool GetNextAlignmentGroup(vector<BamAlignment> &alg);
	bool GetNextAlignmentGroup(BamAlignment &alignment, vector<XATag> &tags);
        const RefVector & GetReferences() const;
        int GetReferenceCount() const;

private:
	void parseXATag(const string &str, vector<XATag> &tags);
	XATag parseXAEntry(const string &str);

private:
	BamReader reader;
	BamAlignment buffer;
};
#endif
