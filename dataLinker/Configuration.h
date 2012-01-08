/*
 * dataLinker : creates abstract contig links from the available information sources.
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

#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H
#include <string>
#include <vector>
#include "AlignerConfiguration.h"

using namespace std;

class PairedInput
{
public:
	PairedInput(const string &leftFileName, const string &rightFileName, double mean, double std, bool isIllumina, double weight = 1, int mapQ = 0, int minReadLength = 0, int maxEditDistance = 10000) : LeftFileName(leftFileName), RightFileName(rightFileName), Mean(mean), Std(std), IsIllumina(isIllumina), Weight(weight), MapQ(mapQ), MinReadLength(minReadLength), MaxEditDistance(maxEditDistance) {};

public:
	string LeftFileName;
	string RightFileName;
	double Mean;
	double Std;
	bool IsIllumina;
	double Weight;
	int MapQ;
	int MinReadLength;
	int MaxEditDistance;
};

class SequenceInput
{
public:
	SequenceInput(const string &fileName, double std, double weight = 1, int minAlignmentLength = 0)
                : FileName(fileName), Std(std), Weight(weight), MinAlignmentLength(minAlignmentLength) {};

public:
	string FileName;
	double Std;
	double Weight;
	int MinAlignmentLength;
};

class Configuration
{
public:
	Configuration();

public:
	bool ProcessCommandLine(int argc, char *argv[]);

public:
	bool Success;
	string InputFileName;
	string OutputFileName;
	string ReadCoverageFileName;
        int MaximumLinkHits;
	double NoOverlapDeviation;
	BWAConfiguration BWAConfig;
	NovoAlignConfiguration NovoAlignConfig;
	SAMToolsConfiguration SAMToolsConfig;
        MummerTilerConfiguration MummerTilerConfig;
	vector<PairedInput> PairedReadInputs;
        vector<SequenceInput> SequenceInputs;
	string LastError;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
