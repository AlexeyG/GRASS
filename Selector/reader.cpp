#include "globals.h"
#include "reader.h"
#include <stdexcept>
#include <iostream>
#include <cstring>

using namespace std;

Reader::Reader()
{
    fin = NULL;
    line = new char[MaxLine];
    buf = new char[MaxLine];
	num_reads = -1;
}

bool Reader::Open(const string &filename, const string &mode)
{
	if (fin != NULL)
		return false;
	fin = fopen(filename.c_str(), mode.c_str());
    if (fin == NULL)
		return false;
	return true;
}

bool Reader::Close()
{
	if (fin != NULL)
	{
        fclose(fin);
		fin = NULL;
		num_reads = -1;
		return true;
	}
	return false;
}

Reader::~Reader()
{
    Close();
    delete [] line;
    delete [] buf;
}

bool FastAReader::Read(string &seq, string &comment)
{
    if (fgets(line, MaxLine, fin) == NULL)
        return false;

    if (line[0] != '>')
        throw runtime_error("FastAReader : comment is not in correct format.");

    line[strlen(line) - 1] = '\0';
    comment = line + 1;
	seq.clear();

    int offset = 0;
    while (fgets(line, MaxLine, fin) != NULL)
    {
        if (line[0] == '>')
        {
            fseek(fin, -strlen(line), SEEK_CUR);
            break;
        }

        int len = strlen(line);
        while (len > 0 && isspace(line[len - 1]))
            --len;
        line[len] = '\0';

        if (len + offset + 1 > MaxLine)
        {
            seq += buf;
            offset = 0;
        }

        strcpy(buf + offset, line);
        offset += len;
    }

    if (offset > 0)
    {
        seq += buf;
        offset = 0;
    }

    for (int i = 0; i < (int)seq.length(); ++i)
        seq[i] = toupper(seq[i]);

    return true;
}

bool FastAReader::Read(FastASequence &seq)
{
	return Read(seq.Nucleotides, seq.Comment);
}

long long FastAReader::Read(vector<FastASequence> &sequences)
{
    sequences.resize(NumReads());
    string seq;
    string comment;
    for (unsigned i = 0; i < sequences.size(); ++i)
    {
        Read(seq, comment);
		sequences[i] = FastASequence(seq, comment);
    }

    return sequences.size();
}

long long FastAReader::NumReads()
{
    long long index = ftell(fin);
    fseek(fin, 0, SEEK_SET);

	if (num_reads >= 0)
		return num_reads;

    long long num = 0;
    while (fgets(line, MaxLine, fin) != NULL)
    {
        if (line[0] == '>')
            ++num;
    }

    fseek(fin, index, SEEK_SET);
    return num_reads = num;
}

bool FastQReader::Read(string &seq, string &comment)
{
	string quality;
	return Read(seq, comment, quality);
}

bool FastQReader::Read(string &seq, string &comment, string &quality)
{
    if (fgets(line, MaxLine, fin) == NULL)
        return false;

    if (line[0] != '@')
		throw runtime_error("FastQReader : comment is not in correct format.");

    line[strlen(line) - 1] = '\0';
    comment = line + 1;
    seq.clear();

    int offset = 0;
    while (fgets(line, MaxLine, fin) != NULL)
    {
        if (line[0] == '+')
        {
			// read the comment
            fgets(line, MaxLine, fin);
            line[strlen(line) - 1] = '\0';
			// read single line (must be shorter than MaxLine) of quality scores
			int len = strlen(line);
			while (len > 0 && isspace(line[len - 1]))
				--len;
			line[len] = '\0';
			quality += line;

            break;
        }

        int len = strlen(line);
        while (len > 0 && isspace(line[len - 1]))
            --len;
        line[len] = '\0';

        if (len + offset + 1 > MaxLine)
        {
            seq += buf;
            offset = 0;
        }

        strcpy(buf + offset, line);
        offset += len;
    }

    if (offset > 0)
    {
        seq += buf;
        offset = 0;
    }

    for (int i = 0; i < (int)seq.length(); ++i)
        seq[i] = toupper(seq[i]);

    return true;
}

bool FastQReader::Read(FastASequence &seq)
{
	return Read(seq.Nucleotides, seq.Comment);
}

bool FastQReader::Read(FastQSequence &seq)
{
	return Read(seq.Nucleotides, seq.Comment, seq.Quality);
}

long long FastQReader::Read(vector<FastASequence> &sequences)
{
    sequences.resize(NumReads());
    string seq;
    string comment;
    for (unsigned i = 0; i < sequences.size(); ++i)
    {
        Read(seq, comment);
		sequences[i] = FastASequence(seq, comment);
    }

    return sequences.size();
}

long long FastQReader::Read(vector<FastQSequence> &sequences)
{
    sequences.resize(NumReads());
    string seq;
    string comment;
	string quality;
    for (unsigned i = 0; i < sequences.size(); ++i)
    {
        Read(seq, comment, quality);
		sequences[i] = FastQSequence(seq, comment, quality);
    }

    return sequences.size();
}

long long FastQReader::NumReads()
{
    long long index = ftell(fin);
    fseek(fin, 0, SEEK_SET);

	if (num_reads >= 0)
		return num_reads;

    long long num = 0;
    while (fgets(line, MaxLine, fin) != NULL)
    {
        if (line[0] == '@')
            ++num;
    }

    fseek(fin, index, SEEK_SET);
    return num_reads = num;
}
