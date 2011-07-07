#ifndef _XATAG_H
#define _XATAG_H
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include <string>

using namespace std;
using namespace BamTools;

class XATag
{
public:
	XATag(int position, int refID, bool isReverseStrand) : Position(position), RefID(refID), IsReverseStrand(isReverseStrand) {};
	XATag(const BamAlignment &alg);
	XATag(const string &str, const BamReader &reader);

public:
	int Position;
	int RefID;
	bool IsReverseStrand;
};
#endif
