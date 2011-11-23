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

int main(int argc, char *argv[])
{
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
