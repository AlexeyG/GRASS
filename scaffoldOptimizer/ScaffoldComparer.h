#ifndef _SCAFFOLDCOMPARER_H
#define _SCAFFOLDCOMPARER_H
#include "ScaffoldExtractor.h"

class ScaffoldComparer
{
public:
	static int Compare(const Scaffold &a, const Scaffold &b);
	static int Compare(const vector<Scaffold> &a, const vector<Scaffold> &b);
	static int OrientationDistance(const Scaffold &a, const Scaffold &b);
	static int OrientationDistance(const vector<Scaffold> &a, const vector<Scaffold> &b);

private:
	static int compareOriented(const Scaffold &a, const Scaffold &b);
};

#endif
