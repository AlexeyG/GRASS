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
