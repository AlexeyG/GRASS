/* 
 * File:   ReadCoverageRepeatDetecter.h
 * Author: alexeyg
 *
 * Created on November 30, 2011, 9:28 AM
 */

#ifndef _READCOVERAGEREPEATDETECTER_H
#define	_READCOVERAGEREPEATDETECTER_H

#include <vector>
#include "ReadCoverage.h"
#include "DataStore.h"

using namespace std;

class ReadCoverageRepeatDetecter
{   
public:
    static vector<int> Detect(double expectedCoverage, const ReadCoverage &coverage, const DataStore &store, double uniqunessCutoff);
};

#endif	/* _READCOVERAGEREPEATDETECTER_H */
