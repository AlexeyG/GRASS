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

/* 
 * File:   ReadCoverageWriter.cpp
 * Author: alexeyg
 * 
 * Created on November 29, 2011, 4:17 PM
 */

#include "ReadCoverageWriter.h"
#include <cstdio>

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
    fprintf(out, "%.10lf\t%lld\t%i\t%i\n", coverage.AverageReadLength, coverage.TotalReadLength, coverage.TotalReadCount, contigCount);
    for (int i = 0; i < contigCount; i++)
    {
        int readCount = (int)coverage.ReadLocations[i].size();
        fprintf(out, "%i\t%i\n", i, readCount);
        for (int j = 0; j < readCount; j++)
            fprintf(out, (j < readCount - 1 ? "%i\t" : "%i\n"), coverage.ReadLocations[i][j]);
    }
    return true;
}
