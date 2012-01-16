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

#include "Writer.h"
#include "Globals.h"
#include "MinMax.h"
#include <string>
#include <cstring>

using namespace std;

Writer::Writer()
{
    fout = NULL;
	buf = new char[MaxLine];
}

Writer::~Writer()
{
    Close();
	delete [] buf;
}

bool Writer::Open(const string &filename, const string &mode)
{
	if (fout != NULL)
		return false;

	fout = fopen(filename.c_str(), mode.c_str());
    if (fout == NULL)
		return false;
	return true;
}

bool Writer::Close()
{
	if (fout != NULL)
	{
        fclose(fout);
		fout = NULL;
		return true;
	}
	return false;
}

void Writer::splitPrint(const string &seq, int num)
{
	int len = seq.length();
	if (num > MaxLine - 1)
		num = MaxLine - 1;
	for (int i = 0; i < len; i += num)
	{
		int left = min(num, len - i);
		memcpy(buf, seq.c_str() + i, left);
		buf[left] = '\0';
		fprintf(fout, "%s\n", buf);
	}
}

bool FastAWriter::Write(const string &seq, const string &comment)
{
	if (fout == NULL)
		return false;

    fprintf(fout, ">%s\n", comment.c_str());
    //fprintf(fout, "%s\n", seq.c_str());
    splitPrint(seq);

    return true;
}

bool FastAWriter::Write(const FastASequence &seq)
{
	return Write(seq.Nucleotides, seq.Comment);
}

bool FastAWriter::Write(const vector<FastASequence> &seq)
{
	int n = seq.size();
	for (int i = 0; i < n; i++)
		if (!Write(seq[i]))
			return false;

	return true;
}

bool FastQWriter::Write(const string &seq, const string &comment)
{
    string quality;
    quality.resize(seq.length());
	fill(quality.begin(), quality.end(), 'a');

    return Write(seq, comment, quality);
}

bool FastQWriter::Write(const string &seq, const string &comment, const string &quality)
{
    if (fout == NULL)
        return false;
    
    fprintf(fout, "@%s\n", comment.c_str());
    fprintf(fout, "%s\n", seq.c_str());
    //splitPrint(seq);
    fprintf(fout, "+%s\n", comment.c_str());
    fprintf(fout, "%s\n", quality.c_str());

    return true;
}

bool FastQWriter::Write(const FastASequence &seq)
{
	return Write(seq.Nucleotides, seq.Comment);
}

bool FastQWriter::Write(const FastQSequence &seq)
{
	return Write(seq.Nucleotides, seq.Comment, seq.Quality);
}

bool FastQWriter::Write(const vector<FastASequence> &seq)
{
	int n = seq.size();
	for (int i = 0; i < n; i++)
		if (!Write(seq[i]))
			return false;
	return true;
}

bool FastQWriter::Write(const vector<FastQSequence> &seq)
{
	int n = seq.size();
	for (int i = 0; i < n; i++)
		if (!Write(seq[i]))
			return false;
	return true;
}
