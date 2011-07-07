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
	//removeBamFiles();
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
			addLinkForTagPair(groupId, *l, leftAlg.Length, *r, rightAlg.Length, input, combinations);
		}
}

void PairedReadConverter::addLinkForTagPair(int groupId, const XATag &l, int lLen, const XATag &r, int rLen, const PairedInput &input, int factor)
{
	int lRefLen = dataStore[l.RefID].Sequence.Nucleotides.length();
	int rRefLen = dataStore[r.RefID].Sequence.Nucleotides.length();
	bool equalOrientation = (l.IsReverseStrand ^ r.IsReverseStrand ? input.IsIllumina : !input.IsIllumina);
	bool forwardOrder = (input.IsIllumina ? !l.IsReverseStrand : l.IsReverseStrand);
	double distance = input.Mean;
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
			distance += l.Position + r.Position;
		else if (!l.IsReverseStrand && r.IsReverseStrand)
			distance += l.Position + rRefLen - r.Position - rLen;
		else if (l.IsReverseStrand && !r.IsReverseStrand)
			distance += lRefLen - l.Position - lLen + r.Position;
		else if (l.IsReverseStrand && r.IsReverseStrand)
			distance += lRefLen - l.Position - lLen + rRefLen - r.Position - rLen;
	}
	//fprintf(stdout, "Test: %c%i %c%i => %s\n", (l.IsReverseStrand ? '-' : '+'), l.RefID, (r.IsReverseStrand ? '-' : '+'), r.RefID, (equalOrientation ? "equal" : "opposite"));
	ContigLink link(l.RefID, r.RefID, distance, input.Std, equalOrientation, forwardOrder, input.Weight / (double)factor);
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
