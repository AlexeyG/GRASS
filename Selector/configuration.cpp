#include "configuration.h"

// Configuration segment class ctor
Segment::Segment(int start, int finish, int chromosome) : Start(start), Finish(finish), Chromosome(chromosome)
{
}

// Comparison operator
bool Segment::operator< (const Segment &other) const
{
	if (Chromosome < other.Chromosome)
		 return true;
	if (Chromosome > other.Chromosome)
		return false;
	if (Start < other.Start)
		return true;
	if (Start > other.Start)
		return false;
	if (Finish > other.Finish)
		return true;
	if (Finish < other.Finish)
		return false;
	return false;
}

// Configuration selection class ctor
ConfigSelect::ConfigSelect(int length, int chromosome, int count) : Length(length), Chromosome(chromosome), Count(count)
{
}

// Configuration input bam class ctor
PairedBam::PairedBam(const string &inputBam1, const string &inputBam2, const string &prefix) : InputBam1(inputBam1), InputBam2(inputBam2), OutputPrefix(prefix)
{
}

// Configuration class defualt ctor
Configuration::Configuration()
{
	InputFastaFileName = "";
	OutputFastaFileName = "";
	PrintChromosomeInfo = false;
	Select.clear();
	Segments.clear();
	PairedAlignment.clear();
}
