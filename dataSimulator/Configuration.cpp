#include "Configuration.h"

Configuration::Configuration()
{
	Success = false;
	AllowContigOverlap = false;
	FlipOrientation = true;
	ShuffleContigs = true;
	Limit = 1;
	Splits.clear();
	InputFileName = "";
	OutputFileName = "out.fasta";
}
