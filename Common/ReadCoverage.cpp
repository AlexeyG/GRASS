#include "ReadCoverage.h"

ReadCoverage::ReadCoverage(int contigCount = 0)
        : ReadLocations(0)
{
    TotalReadCount = 0;
    AverageReadLength = 0;
    totalReadLength = 0;
}

int ReadCoverage::GetContigCount() const
{
    return ReadLocations.size();
}

void ReadCoverage::SetContigCount(int size)
{
    ReadLocations.resize(size, vector<int>());
}

void ReadCoverage::AddLocation(int id, int location)
{
    ReadLocations[id].push_back(location);
}

void ReadCoverage::UpdateAverage(int readLength)
{
    totalReadLength += readLength;
    TotalReadCount++;
    AverageReadLength = (double)totalReadLength / (double)TotalReadCount;
}
