#ifndef _CONTIGINFORMATION_H
#define _CONTIGINFORMATION_H
#include <string>

using namespace std;

class ContigInformation
{
public:
	ContigInformation(int position = 0, bool reverseOrientation = false, int sourceContig = 0) : Position(position), ReverseOrientation(reverseOrientation), SourceContig(sourceContig) {};

public:
	string FormatName() const;

public:
	int Position;
	bool ReverseOrientation;
	int SourceContig;
};
#endif
