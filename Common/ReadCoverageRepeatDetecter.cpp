#include "ReadCoverageRepeatDetecter.h"
#include <algorithm>
#include <cmath>

#include <cstdio>

using namespace std;

vector<int> ReadCoverageRepeatDetecter::Detect(double expectedCoverage, const ReadCoverage &coverage, const DataStore &store, double uniquenessCutoff)
{
    double expectedStarts = expectedCoverage / coverage.AverageReadLength;
    printf("rho: %.5lf\n", expectedStarts);
    vector<int> repeatContigs;
    int nContigs = coverage.ReadLocations.size();
    for (int i = 0; i < nContigs; i++)
    {
        double observedMean = 0.0;
        int contigLength = store[i].Sequence.Nucleotides.length();
        vector<int> readPositions(coverage.ReadLocations[i]);
        sort(readPositions.begin(), readPositions.end());
        vector<int>::const_iterator pos = readPositions.begin();
        while (pos != readPositions.end())
        {
            vector<int>::const_iterator start = pos++;
            int count = 1;
            while (pos != readPositions.end() && *pos == *start)
                count++, pos++;
            if (*start >= 0)
            {
                if (*start <= contigLength - coverage.AverageReadLength)
                    observedMean += count;
            }
        }
        contigLength -= (int)coverage.AverageReadLength;
        if (contigLength <= 0)
            continue;
        observedMean /= (double)contigLength;
        double logRatio = log(2.0) / 2.0 + contigLength * (expectedStarts * expectedStarts - observedMean * observedMean / 2.0) / (2.0 * expectedStarts);
        if (logRatio < uniquenessCutoff)
            repeatContigs.push_back(i);
        if (logRatio < uniquenessCutoff)
            printf("n = %7i, x = %10.2lf, %10.2lf <  %10.2lf -> 0\n", contigLength, observedMean, logRatio, uniquenessCutoff);
        else
            printf("n = %7i, x = %10.2lf, %10.2lf >= %10.2lf -> 1\n", contigLength, observedMean, logRatio, uniquenessCutoff);
    }
    
    return repeatContigs;
}
