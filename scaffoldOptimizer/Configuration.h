/*
 * scaffoldOptimizer : solves the MIQP optimization and produces linear scaffold
 * sequences.
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
#include <sstream>
#include "SolverConfiguration.h"
#include "OverlapperConfiguration.h"

using namespace std;

class Configuration
{
public:
	Configuration();

public:
	bool ProcessCommandLine(int argc, char *argv[]);

public:
	bool Success;
	string LastError;
	bool RemoveAmbiguous;
	bool Sort;
	bool Bundle;
	bool BundlePerGroup;
	bool BundleAmbiguous;
	double BundleDistance;
	double Erosion;
        double ExpectedCoverage;
        double UniquenessFCutoff;
	bool PrintMatrix;
        OverlapperConfiguration OverlapperOptions;
	SolverConfiguration Options;
	string InputFileName;
        string ReadCoverageFileName;
	string OutputFileName;
	string SolutionOutputFileName;

private:
	void printHelpMessage(stringstream &serr);
};

#endif
