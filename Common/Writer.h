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


#ifndef _WRITER_H
#define _WRITER_H

#include <string>
#include <vector>
#include "Sequence.h"

using namespace std;

class Writer
{
public:
    Writer();
    virtual ~Writer();

	bool Open(const string &filename, const string &mode = "wb");
	bool Close();

    virtual bool Write(const string &seq, const string &comment) = 0;
    virtual bool Write(const FastASequence &seq) = 0;
	virtual bool Write(const vector<FastASequence> &seq) = 0;

protected:
	void splitPrint(const string &seq, int num = 80);

protected:
    FILE *fout;
	char *buf;
};

class FastAWriter: public Writer
{
public:
    FastAWriter() : Writer() {}
    ~FastAWriter() {}

    bool Write(const string &seq, const string &comment);
    bool Write(const FastASequence &seq);
	bool Write(const vector<FastASequence> &seq);
};

class FastQWriter: public Writer
{
public:
    FastQWriter() : Writer() {}
    ~FastQWriter() {}

    bool Write(const string &seq, const string &comment);
	bool Write(const string &seq, const string &comment, const string &quality);
    bool Write(const FastASequence &seq);
	bool Write(const vector<FastASequence> &seq);
	bool Write(const FastQSequence &seq);
	bool Write(const vector<FastQSequence> &seq);
};
#endif
