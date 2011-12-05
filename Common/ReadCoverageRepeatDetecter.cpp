#include "ReadCoverageRepeatDetecter.h"
#include <algorithm>
#include <cmath>

#include <cstdio>

using namespace std;

vector<int> ReadCoverageRepeatDetecter::Detect(double expectedCoverage, const ReadCoverage &coverage, const DataStore &store, double uniquenessCutoff)
{
    double expectedStarts = expectedCoverage / coverage.AverageReadLength;
    printf("Expecting coverage: %.5lf\n", expectedStarts);
    vector<int> repeatContigs;
    int nContigs = coverage.ReadLocations.size();
    for (int i = 0; i < nContigs; i++)
    {
        double observedMean = 0.0;
        int observedGroups = 0;
        int contigLength = store[i].Sequence.Nucleotides.length();
        vector<int> readPositions(coverage.ReadLocations[i]);
        sort(readPositions.begin(), readPositions.end());
        vector<int>::const_iterator pos = readPositions.begin();
        while (pos != readPositions.end())
        {
            vector<int>::const_iterator start = pos++;
            int count = 1;
            if (*start >= 0) observedGroups++;
            while (pos != readPositions.end() && *pos == *start)
                count++, pos++;
            if (*start >= 0) observedMean += count;
        }
        //observedMean /= (double)observedGroups;
        observedMean /= (double)contigLength;
        double logRatio = log(2.0) / 2.0 + contigLength * (expectedStarts * expectedStarts - observedMean * observedMean / 2.0) / (2.0 * expectedStarts);
        if (logRatio < uniquenessCutoff)
            repeatContigs.push_back(i);
        printf("Observed coverage: %.5lf -> %.5lf -> %i\n", observedMean, logRatio, logRatio >= uniquenessCutoff);
    }
    
    return repeatContigs;
}
