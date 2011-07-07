#ifndef _CONTIGINFORMATION_H
#define _CONTIGINFORMATION_H
#include <string>

using namespace std;

class ContigInformation
{
public:
	ContigInformation(int position = 0, bool reverseOrientation = false) : Position(position), ReverseOrientation(reverseOrientation) {};

public:
	string FormatName() const;

public:
	int Position;
	bool ReverseOrientation;
};
#endif
