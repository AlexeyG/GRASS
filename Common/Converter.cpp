#include "Converter.h"
#include "Globals.h"
#include "Helpers.h"
#include <string>

#include <iostream>

using namespace std;

bool Converter::Convert(const string &outputFileName, bool removeAfterConversion)
{
	if (InputFileName.length() == 0)
		return false;
	char str[MaxLine];
	OutputFileName = (outputFileName.length() > 0 ? outputFileName : Helpers::TempFile(Configuration.TmpPath));
	sprintf(str, Configuration.ConvertCommand.c_str(), OutputFileName.c_str(), InputFileName.c_str());
	bool result = Helpers::Execute(str);
	if (!result)
	{
		Helpers::RemoveFile(OutputFileName);
		OutputFileName.clear();
	}
	if (removeAfterConversion)
		Helpers::RemoveFile(InputFileName);

	return result;
}

Converter::~Converter()
{
	if (RemoveOutput && OutputFileName.length() > 0)
		Helpers::RemoveFile(OutputFileName);
}
