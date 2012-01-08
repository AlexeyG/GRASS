#ifndef _ALIGNER_H
#define _ALIGNER_H
#include "AlignerConfiguration.h"
#include <string>
#include <vector>

using namespace std;

class Aligner
{
public:
	Aligner(const string &referenceFile, const string &queryFile)
		: ReferenceFileName(referenceFile), QueryFileName(queryFile), RemoveOutput(true) {};
	virtual ~Aligner();

public:
	virtual bool Align(const string &outFile = (char *)"") = 0;

public:
	const string ReferenceFileName;
	const string QueryFileName;
	string OutputFileName;
	bool RemoveOutput;
};

class BWAAligner : public Aligner
{
public:
	BWAAligner(const string &referenceFile, const string &queryFile, const BWAConfiguration &config);

public:
	bool Align(const string &outFile = (char *)"");

public:
	BWAConfiguration Configuration;

private:
	void removeIndexFiles(const string &prefix);

private:
	static vector<string> IndexFileExtensions;
};

class NovoAlignAligner : public Aligner
{
public:
	NovoAlignAligner(const string &referenceFile, const string &queryFile, const NovoAlignConfiguration &config)
		: Aligner(referenceFile, queryFile), Configuration(config) {};

public:
	bool Align(const string &outFile = (char *)"");

public:
	NovoAlignConfiguration Configuration;
};

class MummerAligner : public Aligner
{
public:
    MummerAligner(const string &referenceFile, const string &queryFile, const MummerConfiguration &config)
            : Aligner(referenceFile, queryFile), Configuration(config) {};
          
public:
    bool Align(const string &outFile = (char *)"");
        
public:
    MummerConfiguration Configuration;
};

class MummerTiler : public Aligner
{
public:
    MummerTiler(const string &referenceFile, const string &queryFile, const MummerTilerConfiguration &config)
            : Aligner(referenceFile, queryFile), Configuration(config) {};
          
public:
    bool Align(const string &outFile = (char *)"");
        
public:
    MummerTilerConfiguration Configuration;
};

#endif
