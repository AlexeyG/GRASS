#ifndef _CONTIGINFO_H
#define _CONTIGINFO_H

#include "Sequence.h"

class ContigInfo
{
public:
	ContigInfo(const FastASequence &contig);

public:
	bool operator< (const ContigInfo &other) const;

public:
	FastASequence Contig;
	int Position;
	int Length;

private:
	int determineContigPosition();
};
#endif
