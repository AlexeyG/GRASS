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
