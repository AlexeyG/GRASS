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
 * File:   ReadCoverage.h
 * Author: alexeyg
 *
 * Created on November 29, 2011, 3:41 PM
 */

#ifndef _READCOVERAGE_H
#define	_READCOVERAGE_H

#include <vector>

using namespace std;

class ReadCoverage
{
public:
    ReadCoverage(int contigCount = 0);
    
public:
    int GetContigCount() const;
    void SetContigCount(int count);
    void AddLocation(int id, int location);
    void UpdateAverage(int readLength);
    void SetAverageReadLength(long long totalReadLength, int totalReadCount);
    
public:
    vector< vector<int> > ReadLocations;
    double AverageReadLength;
    int TotalReadCount;
    long long TotalReadLength;
};

#endif	/* READCOVERAGE_H */
