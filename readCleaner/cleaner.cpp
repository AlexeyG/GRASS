/*
 * readCleaner : a tool used in debugging. Removes paired reads originating from
 * gaps
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

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include "Configuration.h"
#include "PairedReadProcessor.h"

using namespace std;

Configuration config;

bool processPairs(const Configuration &config, const vector<PairedInput> &paired)
{
	PairedReadProcessor processor;
	int n = paired.size();
	for (int i = 0; i < n; i++)
	{
		PairedInput p = paired[i];
		cerr << "   [i] Processing " << (p.IsIllumina ? "Illumina" : "454") << " paired reads (" << p.LeftFileName << ", " << p.RightFileName << ")" << endl; 
		switch (processor.Process(config, paired[i]))
		{
		case PairedReadProcessor::Success:
			cerr << "      [+] Successfully processed paired reads." << endl;
			break;
		case PairedReadProcessor::FailedLeftAlignment:
			cerr << "      [-] Unable to align left read mates (" << (p.IsIllumina ? "Illumina" : "454") << ")." << endl;
			return false;
		case PairedReadProcessor::FailedRightAlignment:
			cerr << "      [-] Unable to align second read mates (" << (p.IsIllumina ? "Illumina" : "454") << ")." << endl;
			return false;
		case PairedReadProcessor::FailedLeftConversion:
			cerr << "      [-] Unable to convert left alignment SAM to BAM (SAM Tools)." << endl;
			return false;
		case PairedReadProcessor::FailedRightConversion:
			cerr << "      [-] Unable to convert right alignment SAM to BAM (SAM Tools)." << endl;
			return false;
		case PairedReadProcessor::FailedIO:
			cerr << "      [-] Unable to create open alignments or output filtered reads (IO failure)." << endl;
			return false;
		}
	}
	return true;
}

void banner()
{
    cerr << "This program comes with ABSOLUTELY NO WARRANTY; see LICENSE for details." << endl;
    cerr << "This is free software, and you are welcome to redistribute it" << endl;
    cerr << "under certain conditions; see LICENSE for details." << endl;
    cerr << endl;
}

int main(int argc, char *argv[])
{
    banner();
	srand((unsigned int)time(NULL));
	if (config.ProcessCommandLine(argc, argv))
	{
		cerr << "[i] Processing paired reads." << endl;
		if (!processPairs(config, config.PairedReadInputs))
			return -1;
		return 0;
	}
	cerr << config.LastError;
	return -1;
}
