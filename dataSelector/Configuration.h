/*
 * dataSelector : a support tool used in debugging to create simulated scaffolding
 * problems given a reference sequence and paired reads.
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

class PairedSimulation
{
public:
	PairedSimulation(double readLengthMean, double readLengthStd, double insertSizeMean, double insertSizeStd, double depth, bool isIllumina, const string &prefix);

public:
	double ReadLengthMean;
	double ReadLengthStd;
	double InsertSizeMean;
	double InsertSizeStd;
	double Depth;
	bool IsIllumina;
	string OutputPrefix;
};

class Configuration
{
public:
	Configuration();

public:
	bool ProcessCommandLine(int argc, char *argv[]);
	bool CheckConfig() const;

public:
	bool Success;
	string LastError;
	string InputFastaFileName;
	string OutputFastaFileName;
	bool PrintChromosomeInfo;
	vector<Segment> Segments;
	vector<ConfigSelect> Select;
	vector<PairedBam> PairedAlignment;
	vector<PairedSimulation> PairedReadSimulation;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
