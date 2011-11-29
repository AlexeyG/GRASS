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
