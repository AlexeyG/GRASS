#ifndef _SCAFFOLDCONVERTER_H
#define _SCAFFOLDCONVERTER_H

#include "ScaffoldExtractor.h"
#include "DataStore.h"
#include "Sequence.h"
#include "OverlapperConfiguration.h"
#include <vector>

using namespace std;

class ScaffoldConverter
{
public:
	static FastASequence ToFasta(const DataStore &store, const Scaffold &scaffold, const OverlapperConfiguration &config);
	static vector<FastASequence> ToFasta(const DataStore &store, const vector<Scaffold> &scaffold, const OverlapperConfiguration &config);
};

#endif
