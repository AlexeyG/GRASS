#ifndef _CONVERTER_H
#define _CONVERTER_H

#include "Configuration.h"
#include <string>

using namespace std;

class Converter
{
public:
	Converter(const string &inputFileName, const SAMToolsConfiguration &config)
		: InputFileName(inputFileName), Configuration(config), RemoveOutput(true) {};
	virtual ~Converter();

public:
	bool Convert(const string &outputFileName = (char *)"" ,bool removeAfterConversion = true);

public:
	const string InputFileName;
	string OutputFileName;
	SAMToolsConfiguration Configuration;
	bool RemoveOutput;
};
#endif
