#include "PairedReadConverter.h"
#include "Aligner.h"
#include "Converter.h"
#include "Helpers.h"
#include "AlignmentReader.h"
#include <sstream>

PairedReadConverter::PairedReadConverter(DataStore &store)
	: dataStore(store) 
{
}

bool PairedReadConverter::IsCorrectRelativeOrientation(const XATag &l, const XATag &r, bool isIllumina)
{
	return (isIllumina ? l.IsReverseStrand != r.IsReverseStrand : l.IsReverseStrand == r.IsReverseStrand);
}

PairedReadConverter::PairedReadConverterResult PairedReadConverter::Process(const Configuration &config, const PairedInput &input)
{
	PairedReadConverterResult result = alignAndConvert(config, input);
	if (result == Success)
	{
		string groupName = (input.IsIllumina ? "Ilumina paired read alignment" : "454 paired read alignment");
		stringstream groupDescription;
		groupDescription << (input.IsIllumina ? "Illumina" : "454") << " paired reads: " << input.LeftFileName << " & " << input.RightFileName << " with " << input.Mean << " +/- " << input.Std << " of weight " << input.Weight;
		int groupId = dataStore.AddGroup(LinkGroup(groupName, groupDescription.str()));
		result = createLinksFromAlignment(groupId, config.MaximumLinkHits, input);
	}
	removeBamFiles();
	return result;
}

PairedReadConverter::PairedReadConverterResult PairedReadConverter::alignAndConvert(const Configuration &config, const PairedInput &input)
{
	Aligner *leftAlignment, *rightAlignment;
	PairedReadConverterResult result = Success;
	bool success = true;

	leftAlignment = (input.IsIllumina ? (Aligner *)new BWAAligner(config.InputFileName, input.LeftFileName, config.BWAConfig) : new NovoAlignAligner(config.InputFileName, input.LeftFileName, config.NovoAlignConfig));
	rightAlignment = (input.IsIllumina ? (Aligner *)new BWAAligner(config.InputFileName, input.RightFileName, config.BWAConfig) : new NovoAlignAligner(config.InputFileName, input.RightFileName, config.NovoAlignConfig));

	success = success && leftAlignment->Align();
	if (!success)
		result = FailedLeftAlignment;
	else
	{
		success = success && rightAlignment->Align();
		if (!success)
			result = FailedRightAlignment;
	}

	Converter leftConverter(leftAlignment->OutputFileName, config.SAMToolsConfig);
	Converter rightConverter(rightAlignment->OutputFileName, config.SAMToolsConfig);
	if (success)
	{
		success = success && leftConverter.Convert();
		if (!success)
			result = FailedLeftConversion;
		else
		{
			success = success && rightConverter.Convert();
			if (!success)
				result = FailedRightConversion;
		}
	}

	if (success)
	{
		leftConverter.RemoveOutput = rightConverter.RemoveOutput = false;
		leftBamFileName = leftConverter.OutputFileName;
		rightBamFileName = rightConverter.OutputFileName;
	}

	delete leftAlignment;
	delete rightAlignment;

	return result;
}

PairedReadConverter::PairedReadConverterResult PairedReadConverter::createLinksFromAlignment(int groupId, int maxHits, const PairedInput &input)
{
	PairedReadConverterResult result = Success;
	AlignmentReader leftReader, rightReader;
	if (!leftReader.Open(leftBamFileName) || !rightReader.Open(rightBamFileName))
		result = FailedLinkCreation;
	vector<XATag> leftTags, rightTags;
	BamAlignment leftAlignment, rightAlignment;
	while (leftReader.GetNextAlignmentGroup(leftAlignment, leftTags) && rightReader.GetNextAlignmentGroup(rightAlignment, rightTags))
		createLinksForPair(groupId, leftAlignment, leftTags, rightAlignment, rightTags, input, maxHits);
	leftReader.Close();
	rightReader.Close();
	return result;
}

void PairedReadConverter::createLinksForPair(int groupId, const BamAlignment &leftAlg, const vector<XATag> &leftTags, const BamAlignment &rightAlg, const vector<XATag> &rightTags, const PairedInput &input, int maxHits)
{
	int combinations = leftTags.size() * rightTags.size();
	if (combinations > maxHits || combinations == 0)
		return;
	for (vector<XATag>::const_iterator l = leftTags.begin(); l != leftTags.end(); l++)
		for (vector<XATag>::const_iterator r = rightTags.begin(); r != rightTags.end(); r++)
		{
			if (l->RefID == r->RefID)
				continue;
			addLinkForTagPair(groupId, *l, leftAlg, *r, rightAlg, input, combinations);
		}
}

void PairedReadConverter::addLinkForTagPair(int groupId, const XATag &l, const BamAlignment &leftAlg, const XATag &r, const BamAlignment &rightAlg, const PairedInput &input, int factor)
{
	int lRefLen = dataStore[l.RefID].Sequence.Nucleotides.length();
	int rRefLen = dataStore[r.RefID].Sequence.Nucleotides.length();
	int lLen = leftAlg.Length;
	int rLen = rightAlg.Length;
	bool equalOrientation = (l.IsReverseStrand ^ r.IsReverseStrand ? input.IsIllumina : !input.IsIllumina);
	bool forwardOrder = l.IsReverseStrand ^ input.IsIllumina;
	double distance = input.Mean;

	if (l.Position + lLen >= lRefLen || r.Position + rLen >= rRefLen)
		return;
	if (leftAlg.MapQuality < input.MapQ || rightAlg.MapQuality < input.MapQ)
		return;
	if (lLen < input.MinReadLength || rLen < input.MinReadLength)
		return;

	/*if (inList(leftAlg.Name) || inList(rightAlg.Name))
	{
		unsigned int distA = leftAlg.GetEditDistance(distA) ? distA : 0;
		unsigned int distB = leftAlg.GetEditDistance(distB) ? distB : 0;
		printf("Left:  %s pos = %4i end = %4i of [%i %i) => %i quality = %i dist = %i\n", leftAlg.Name.c_str(), leftAlg.Position, leftAlg.GetEndPosition(), getStart(dataStore[l.RefID].Sequence.Name(), lRefLen), getEnd(dataStore[l.RefID].Sequence.Name(), lRefLen), getReadPosition(leftAlg), leftAlg.MapQuality, (int)distA);
		printf("Right: %s pos = %4i end = %4i of [%i %i) => %i quality = %i dist = %i\n", rightAlg.Name.c_str(), rightAlg.Position, rightAlg.GetEndPosition(), getStart(dataStore[r.RefID].Sequence.Name(), rRefLen), getEnd(dataStore[r.RefID].Sequence.Name(), rRefLen), getReadPosition(rightAlg), rightAlg.MapQuality, (int)distB);
		printf("\n");
	}
	else
		printf("MAPQ: %i %i\n", leftAlg.MapQuality, rightAlg.MapQuality);*/

	if (!input.IsIllumina)
	{
		if (!l.IsReverseStrand && !r.IsReverseStrand)
			distance += lRefLen - l.Position - lLen + r.Position;
		else if (!l.IsReverseStrand && r.IsReverseStrand)
			distance += lRefLen - l.Position - lLen + rRefLen - r.Position - rLen;
		else if (l.IsReverseStrand && !r.IsReverseStrand)
			distance += l.Position + r.Position;
		else if (l.IsReverseStrand && r.IsReverseStrand)
			distance += l.Position + rRefLen - r.Position - rLen;
	}
	else
	{
		if (!l.IsReverseStrand && !r.IsReverseStrand)
		{
			//printf("A\n");
			distance += l.Position + r.Position;
		}
		else if (!l.IsReverseStrand && r.IsReverseStrand)
		{
			//printf("B\n");
			distance += l.Position + rRefLen - r.Position - rLen;
		}
		else if (l.IsReverseStrand && !r.IsReverseStrand)
		{
			//printf("C\n");
			distance += lRefLen - l.Position - lLen + r.Position;
		}
		else if (l.IsReverseStrand && r.IsReverseStrand)
		{
			//printf("D\n");
			distance += lRefLen - l.Position - lLen + rRefLen - r.Position - rLen;
		}
	}
	//distance += lLen + rLen;
	//fprintf(stdout, "Test: %c%i %c%i => %s\n", (l.IsReverseStrand ? '-' : '+'), l.RefID, (r.IsReverseStrand ? '-' : '+'), r.RefID, (equalOrientation ? "equal" : "opposite"));
	string comment = leftAlg.Name + "+" + rightAlg.Name;
	ContigLink link(l.RefID, r.RefID, distance, input.Std, equalOrientation, forwardOrder, input.Weight / (double)factor, comment);
	link.Ambiguous = factor > 1;
	dataStore.AddLink(groupId, link);
}

void PairedReadConverter::removeBamFiles()
{
	if (!leftBamFileName.empty())
		Helpers::RemoveFile(leftBamFileName);
	if (!rightBamFileName.empty())
		Helpers::RemoveFile(rightBamFileName);
	leftBamFileName.clear();
	rightBamFileName.clear();
}

bool PairedReadConverter::inList(const string &name)
{
	string list = "133448.1|74726+133448.2|75182|286248.1|74732+286248.2|75185|246773.1|74733+246773.2|75186|110685.1|74733+110685.2|75185|219041.1|74728+219041.2|75180|197878.1|74727+197878.2|75178|315396.1|74735+315396.2|75185|278283.1|74727+278283.2|75176|113656.1|74736+113656.2|75184|363369.1|74731+363369.2|75174|177996.1|74740+177996.2|75182|19598.1|74735+19598.2|75175|667.1|74742+667.2|75181|197133.1|74731+197133.2|75169|153451.1|74733+153451.2|75164|85752.1|74757+85752.2|75184|139282.1|74739+139282.2|75162|63446.1|74761+63446.2|75179|237589.1|74743+237589.2|75160|197788.1|74766+197788.2|75174";
	return list.find(name) != string::npos;
}

int PairedReadConverter::getStart(const string &name, int len)
{
	int i = 1, strLen = name.length();
	while (i < strLen && name[i] != '|')
		i++;
	int num = 0;
	for (int j = 1; j < i; j++)
		num = num * 10 + name[j] - '0';
	if (name[0] == '-')
		num -= len;
	return num;
}

int PairedReadConverter::getEnd(const string &name, int len)
{
	return getStart(name, len) + len;
}

int PairedReadConverter::getReadPosition(const BamAlignment &alg)
{
	string name = dataStore[alg.RefID].Sequence.Name();
	int len = dataStore[alg.RefID].Sequence.Nucleotides.length();
	int start = getStart(name, len);
	int pos = alg.Position;
	if (name[0] == '-')
		pos = len - pos - alg.Length;
	return start + pos;
}
