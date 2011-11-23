#ifndef _SCAFFOLDCONVERTER_H
#define _SCAFFOLDCONVERTER_H

#include "ScaffoldExtractor.h"
#include "DataStore.h"
#include "Sequence.h"
#include <vector>

using namespace std;

class ScaffoldConverter
{
public:
	static FastASequence ToFasta(const DataStore &store, const Scaffold &scaffold);
	static vector<FastASequence> ToFasta(const DataStore &store, const vector<Scaffold> &scaffold);
};

#endif
