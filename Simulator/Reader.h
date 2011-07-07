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
