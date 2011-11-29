#include "ReadCoverageReader.h"
#include "Helpers.h"

ReadCoverageReader::ReadCoverageReader()
{
}

ReadCoverageReader::~ReadCoverageReader()
{
    if (in.is_open())
        in.close();
}

bool ReadCoverageReader::Open(const string &fileName)
{
    if (in.is_open())
        return false;
    in.open(fileName.c_str(), ios::in);
    return in.is_open();
}

bool ReadCoverageReader::Close()
{
    in.close();
    return true;
}

bool ReadCoverageReader::Read(ReadCoverage &coverage)
{
    int nContigs;
    if (!in.is_open())
        return false;
    coverage = ReadCoverage();
    if (!readHeader(nContigs, coverage))
        return false;
    if (!readContigs(nContigs, coverage))
        return false;
    return true;
}

bool ReadCoverageReader::readHeader(int &nContigs, ReadCoverage &coverage)
{
    string line;
    getline(in, line);
    string averageReadLengthStr = Helpers::NextEntry(line);
    string totalReadLengthStr = Helpers::NextEntry(line);
    string nReadsStr = Helpers::NextEntry(line);
    string nContigsStr = Helpers::NextEntry(line);
    if (averageReadLengthStr.length() == 0 || nReadsStr.length() == 0 || totalReadLengthStr.length() == 0 || nContigsStr.length() == 0)
        return false;
    double averageReadLength = Helpers::GetArgument<double>(averageReadLengthStr);
    long long totalReadLength = Helpers::GetArgument<long long>(totalReadLengthStr);
    int nReads = Helpers::GetArgument<int>(nReadsStr);
    nContigs = Helpers::GetArgument<int>(nContigsStr);
    if (averageReadLength <= 0 || totalReadLength < 0 || nReads <= 0 || nContigs < 0)
        return false;
    
    coverage.SetContigCount(nContigs);
    coverage.SetAverageReadLength(totalReadLength, nReads);
    
    return true;
}

bool ReadCoverageReader::readContigs(int nContigs, ReadCoverage &coverage)
{
    for (int i = 0; i < nContigs; i++)
        if (!readContig(i, coverage))
            return false;
    return true;
}

bool ReadCoverageReader::readContig(int contigID, ReadCoverage &coverage)
{
    string line;
    getline(in, line);
    string contigIdStr = Helpers::NextEntry(line);
    string nPositionsStr = Helpers::NextEntry(line);
    int nPositions = Helpers::GetArgument<int>(line);
    if (nPositions < 0)
        return false;
    getline(in, line);
    for (int j = 0; j < nPositions; j++)
        coverage.AddLocation(contigID, Helpers::GetArgument<int>(Helpers::NextEntry(line)));
}
