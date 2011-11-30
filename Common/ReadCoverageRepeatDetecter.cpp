#include "ReadCoverageRepeatDetecter.h"
#include <algorithm>
#include <cmath>

#include <cstdio>
using namespace std;

vector<int> ReadCoverageRepeatDetecter::Detect(double expectedCoverage, const ReadCoverage &coverage, const DataStore &store, double uniqunessCutoff)
{
    double expectedStarts = expectedCoverage / coverage.AverageReadLength;
    printf("Expected starts: %.5lf\n", expectedStarts);
    vector<int> repeatContigs;
    int nContigs = coverage.ReadLocations.size();
    for (int i = 0; i < nContigs; i++)
    {
        double observedMean = 0.0;
        int observedGroups = 0;
        int contigLength = store[i].Sequence.Nucleotides.length();
        vector<int> readPositions(coverage.ReadLocations[i]);
        printf("Contig %i: positions %i position %i\n", i, coverage.ReadLocations[i].size(), readPositions.size());
        sort(readPositions.begin(), readPositions.end());
        vector<int>::const_iterator pos = readPositions.begin();
        while (pos != readPositions.end())
        {
            observedGroups++;
            int count = 1;
            vector<int>::const_iterator start = pos++;
            while (pos != readPositions.end() && *pos == *start)
                count++, pos++;
            observedMean += count;
            printf("          Group %i count: %i\n", observedGroups - 1, count);
        }
        printf("          Total groups: %i\n", observedGroups);
        observedMean /= (double)observedGroups;
        printf("Contig %i observed mean: %.5lf\n", i, observedMean);
        double logRatio = log(2.0) / 2.0 + contigLength * (expectedStarts * expectedStarts - observedMean * observedMean / 2.0) / (2.0 * expectedStarts);
        printf("Contig %i log ratio: %.5lf\n", i, logRatio);
        if (logRatio < uniqunessCutoff)
            repeatContigs.push_back(i);
    }
    
    return repeatContigs;
}
