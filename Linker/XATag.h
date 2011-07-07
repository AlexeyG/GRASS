#ifndef _XATAG_H
#define _XATAG_H

#include "api/BamAlignment.h"
#include <string>

using namespace std;
using namespace BamTools;

class XATag
{
public:
	XATag(int position, int refID, bool isReverseStrand) : Position(position), RefID(refID), IsReverseStrand(isReverseStrand) {};
	XATag(const BamAlignment &alg);

public:
	int Position;
	int RefID;
	bool IsReverseStrand;
};
#endif
