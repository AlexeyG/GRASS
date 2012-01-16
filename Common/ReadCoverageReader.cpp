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

#include "ReadCoverageReader.h"
#include "Helpers.h"
#include <sstream>

using namespace std;

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
    coverage.SetContigCount(0);
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
        if (!readContig(coverage))
            return false;
    return true;
}

bool ReadCoverageReader::readContig(ReadCoverage &coverage)
{
    string line;
    getline(in, line);
    string contigIdStr = Helpers::NextEntry(line);
    string nPositionsStr = Helpers::NextEntry(line);
    int contigID = Helpers::GetArgument<int>(contigIdStr);
    int nPositions = Helpers::GetArgument<int>(nPositionsStr);
    if (contigID < 0)
        return false;
    if (nPositions < 0)
        return false;
    getline(in, line);
    stringstream ss(line);
    for (int j = 0, pos; j < nPositions; j++)
    {
        ss >> pos;
        coverage.AddLocation(contigID, pos);
    }
    return true;
}
