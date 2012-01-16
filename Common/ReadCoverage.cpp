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

#include "ReadCoverage.h"

ReadCoverage::ReadCoverage(int contigCount)
        : ReadLocations(contigCount)
{
    TotalReadCount = 0;
    AverageReadLength = 0;
    TotalReadLength = 0;
}

int ReadCoverage::GetContigCount() const
{
    return ReadLocations.size();
}

void ReadCoverage::SetContigCount(int count)
{
    ReadLocations.resize(count, vector<int>());
}

void ReadCoverage::AddLocation(int id, int location)
{
    ReadLocations[id].push_back(location);
}

void ReadCoverage::UpdateAverage(int readLength)
{
    TotalReadLength += readLength;
    TotalReadCount++;
    AverageReadLength = (double)TotalReadLength / (double)TotalReadCount;
}

void ReadCoverage::SetAverageReadLength(long long totalReadLength, int totalReadCount)
{
    TotalReadLength = totalReadLength;
    TotalReadCount = totalReadCount;
    AverageReadLength = (double)totalReadLength / (double)totalReadCount;
}
