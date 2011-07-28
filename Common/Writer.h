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
	void splitPrint(const string &seq, int num = 1000);

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
