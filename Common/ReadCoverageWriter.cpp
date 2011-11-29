/* 
 * File:   ReadCoverageWriter.cpp
 * Author: alexeyg
 * 
 * Created on November 29, 2011, 4:17 PM
 */

#include "ReadCoverageWriter.h"

ReadCoverageWriter::ReadCoverageWriter()
{
    out = NULL;
}

ReadCoverageWriter::~ReadCoverageWriter()
{
    if (out != NULL)
    {
        fclose(out);
        out = NULL;
    }
}

bool ReadCoverageWriter::Open(const string &fileName, const string &mode)
{
    if (out != NULL)
        return false;
    out = fopen(fileName.c_str(), mode.c_str());
    return out != NULL;
}

bool ReadCoverageWriter::Close()
{
    fclose(out);
    out = NULL;
    return true;
}

bool ReadCoverageWriter::Write(const ReadCoverage &coverage)
{
    if (out == NULL)
        return false;
    int contigCount = coverage.GetContigCount();
    fprintf(out, "%.10lf\t%i\t%i\n", coverage.AverageReadLength, coverage.TotalReadCount, contigCount);
    for (int i = 0; i < contigCount; i++)
    {
        int readCount = (int)coverage.ReadLocations[i].size();
        fprintf(out, "%i\t%i\n", i, readCount);
        for (int j = 0; j < readCount; j++)
            fprintf(out, (j < readCount - 1 > "%i\t", "%i\n"), coverage.ReadLocations[i][j]);
    }
    return true;
}
