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
 * Defines a number of helpful functions.
 */

#ifndef _READER_H
#define _READER_H

#include <cstdio>
#include <string>
#include <vector>
#include "Sequence.h"

using namespace std;

class Reader
{
public:
    Reader();
    virtual ~Reader();

    virtual bool Read(FastASequence &seq) = 0;
    virtual bool Read(string &seq, string &comment) = 0;
    virtual long long NumReads() = 0;
    virtual long long Read(vector<FastASequence> &sequences) = 0;
    bool Open(const string &filename, const string &mode = "rb");
    bool Close();
    bool IsOpen() const;
    void Rewind() { fseek(fin, 0, SEEK_SET); }

protected:
    long long num_reads;
    FILE *fin;
    char *line;
    char *buf;
};

class FastAReader: public Reader
{
public:
    FastAReader() : Reader() {}
    virtual ~FastAReader() {}

    bool Read(string &seq, string &comment);
    bool Read(FastASequence &seq);
    long long Read(vector<FastASequence> &sequences);
    long long NumReads();
};

class FastQReader: public Reader
{
public:
    FastQReader() : Reader() {}
    virtual ~FastQReader() {}

    bool Read(string &seq, string &comment);
	bool Read(string &seq, string &comment, string &quality);
	bool Read(FastASequence &seq);
	bool Read(FastQSequence &seq);
	long long Read(vector<FastASequence> &sequences);
	long long Read(vector<FastQSequence> &sequences);
    long long NumReads();
};

#endif
