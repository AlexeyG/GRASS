/*
 * Defines a number of helpful classes and structures.
 */
#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H

#include <string>
#include <vector>

using namespace std;

class Segment
{
public:
	Segment() : Start(-1), Finish(-1), Chromosome(-1) {};
	Segment(int start, int finish, int chromosome);
	bool operator< (const Segment &other) const;
public:
	int Start;
	int Finish;
	int Chromosome;
};

class ConfigSelect
{
public:
	ConfigSelect(int length, int chromosome, int count);

public:
	int Length;
	int Chromosome;
	int Count;
};

class PairedBam
{
public:
	PairedBam(const string &inputBam1, const string &inputBam2, const string &prefix);

public:
	string InputBam1;
	string InputBam2;
	string OutputPrefix;
};

class Configuration
{
public:
	Configuration();

public:
	string InputFastaFileName;
	string OutputFastaFileName;
	bool PrintChromosomeInfo;
	vector<Segment> Segments;
	vector<ConfigSelect> Select;
	vector<PairedBam> PairedAlignment;
};

#endif
