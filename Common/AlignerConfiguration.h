#ifndef _ALIGNERCONFIGURATION_H
#define _ALIGNERCONFIGURATION_H

#include <string>

using namespace std;

class BWAConfiguration
{
public:
	BWAConfiguration();

public:
	int NumberOfThreads;
	int MaximumHits;
	bool ExactMatch;
	string IndexCommand;
	string SuffixArrayCommand;
	string SuffixArrayExactCommand;
	string AlignSingleEndCommand;
	string TmpPath;
};

class NovoAlignConfiguration
{
public:
	NovoAlignConfiguration();

public:
	string IndexCommand;
	string AlignSingleEndCommand;
	string TmpPath;
};

class SAMToolsConfiguration
{
public:
	SAMToolsConfiguration();

public:
	string ConvertCommand;
	string TmpPath;
};

class MummerConfiguration
{
public:
    MummerConfiguration();
    
public:
    string TmpPath;
    string NucmerCommand;
    string DeltaFilterCommand;
    string ShowCoordsCommand;
};

class MummerTilerConfiguration
{
public:
    MummerTilerConfiguration();
    
public:
    string TmpPath;
    string NucmerCommand;
    string ShowTilingCommand;
};

#endif

