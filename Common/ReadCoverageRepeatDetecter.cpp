/*
 * Common : a collection of classes (re)used throughout the scaffolder implementation.
 * Copyright (C) 2011  Alexey Gritsenko
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 * 
 * 
 * Email: a.gritsenko@tudelft.nl
 * Mail: Delft University of Technology
 *       Faculty of Electrical Engineering, Mathematics, and Computer Science
 *       Department of Mediamatics
 *       P.O. Box 5031
 *       2600 GA, Delft, The Netherlands
 */

#include "ReadCoverageRepeatDetecter.h"
#include <algorithm>
#include <cmath>

#include <cstdio>

using namespace std;

vector<int> ReadCoverageRepeatDetecter::Detect(double expectedCoverage, const ReadCoverage &coverage, const DataStore &store, double uniquenessCutoff)
{
    double expectedStarts = expectedCoverage / coverage.AverageReadLength;
    //printf("rho: %.5lf\n", expectedStarts);
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
        /*if (logRatio < uniquenessCutoff)
            printf("n = %7i, x = %10.2lf, %10.2lf <  %10.2lf -> 0\n", contigLength, observedMean, logRatio, uniquenessCutoff);
        else
            printf("n = %7i, x = %10.2lf, %10.2lf >= %10.2lf -> 1\n", contigLength, observedMean, logRatio, uniquenessCutoff);*/
    }
    
    return repeatContigs;
}
