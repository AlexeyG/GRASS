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
 * File:   MummerCoordReader.h
 * Author: alexeyg
 *
 * Created on December 6, 2011, 9:48 AM
 */

#ifndef _MUMMERCOORDREADER_H
#define	_MUMMERCOORDREADER_H

#include <vector>
#include <cstdio>
#include <string>
#include <map>
#include "MummerCoord.h"
#include "Sequence.h"

using namespace std;

class MummerCoordReader
{
public:
    MummerCoordReader(const vector<FastASequence> &references, const vector<FastASequence> &scaffolds);
    ~MummerCoordReader();
    
public:
    bool Open(const string &fileName, const string &mode = "rb");
    bool Close();
    bool IsOpen() const;
    bool Read(MummerCoord &coord);
    long long Read(vector<MummerCoord> &coords);
    long long NumCoords();
    
private:
    FILE *fin;
    char *line;
    long long numCoords;
    
private:
    map<string, int> referenceIds;
    map<string, int> scaffoldIds;
    
private:
    void createMap(map<string, int> &store, const vector<FastASequence> &seq);
};

#endif	/* _MUMMERCOORDREADER_H */
