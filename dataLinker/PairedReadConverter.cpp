#include "PairedReadConverter.h"
#include "Aligner.h"
#include "Converter.h"
#include "Helpers.h"
#include "AlignmentReader.h"
#include <sstream>

using namespace BamTools;

PairedReadConverter::PairedReadConverter(DataStore &store)
	: dataStore(store) 
{
    AverageReadLength = 0;
    TotalReadCount = 0;
    ReadLocations.clear();
    
    totalReadLength = 0;
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
		result = createLinksFromAlignment(groupId, config.MaximumLinkHits, input, config.NoOverlapDeviation);
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
		//cout << input.IsIllumina << " " << input.LeftFileName << " " << leftBamFileName << endl;
		//cout << input.IsIllumina << " " << input.RightFileName << " " << rightBamFileName << endl;
	}

	delete leftAlignment;
	delete rightAlignment;

	return result;
}

PairedReadConverter::PairedReadConverterResult PairedReadConverter::createLinksFromAlignment(int groupId, int maxHits, const PairedInput &input, double noOverlapDeviation)
{
	PairedReadConverterResult result = Success;
	AlignmentReader leftReader, rightReader;
	if (!leftReader.Open(leftBamFileName) || !rightReader.Open(rightBamFileName))
		result = FailedLinkCreation;
	vector<XATag> leftTags, rightTags;
	BamAlignment leftAlignment, rightAlignment;
	while (leftReader.GetNextAlignmentGroup(leftAlignment, leftTags) && rightReader.GetNextAlignmentGroup(rightAlignment, rightTags))
        {
            processCoverage(leftAlignment, leftTags);
            processCoverage(rightAlignment, rightTags);
            createLinksForPair(groupId, leftAlignment, leftTags, rightAlignment, rightTags, input, noOverlapDeviation, maxHits);
            AverageReadLength = (double)totalReadLength / (double)TotalReadCount;
        }
	leftReader.Close();
	rightReader.Close();
	return result;
}

void PairedReadConverter::createLinksForPair(int groupId, const BamAlignment &leftAlg, const vector<XATag> &leftTags, const BamAlignment &rightAlg, const vector<XATag> &rightTags, const PairedInput &input, double noOverlapDeviation, int maxHits)
{
	int combinations = leftTags.size() * rightTags.size();
	/*if (leftAlg.Name == "660005.1|479132" | leftAlg.Name == "660005.2|479603")
	{
		cout << leftAlg.IsMapped() << " " << leftAlg.Length << " " << leftAlg.MapQuality << " " << leftAlg.Name << " " << leftAlg.Position << endl;
		cout << rightAlg.IsMapped() << " " << rightAlg.Length << " " << rightAlg.MapQuality << " " << rightAlg.Name << " " << rightAlg.Position << endl;
	}*/
	if (combinations > maxHits || combinations == 0)
		return;
	for (vector<XATag>::const_iterator l = leftTags.begin(); l != leftTags.end(); l++)
		for (vector<XATag>::const_iterator r = rightTags.begin(); r != rightTags.end(); r++)
		{
			if (l->RefID == r->RefID)
				continue;
			addLinkForTagPair(groupId, *l, leftAlg, *r, rightAlg, input, noOverlapDeviation, combinations);
		}
}

void PairedReadConverter::processCoverage(const BamAlignment &alg, const vector<XATag> &tags)
{
    if (tags.size() != 1)
        return;
    totalReadLength += alg.QueryBases.length();
    TotalReadCount++;
    AverageReadLength = (double) totalReadLength / (double) TotalReadCount;
    
}

void PairedReadConverter::addLinkForTagPair(int groupId, const XATag &l, const BamAlignment &leftAlg, const XATag &r, const BamAlignment &rightAlg, const PairedInput &input, double noOverlapDeviation, int factor)
{
	int lRefLen = dataStore[l.RefID].Sequence.Nucleotides.length();
	int rRefLen = dataStore[r.RefID].Sequence.Nucleotides.length();
	int lLen = leftAlg.Length;
	int rLen = rightAlg.Length;
	bool equalOrientation = (l.IsReverseStrand ^ r.IsReverseStrand ? input.IsIllumina : !input.IsIllumina);
	bool forwardOrder = l.IsReverseStrand ^ input.IsIllumina;
	double distance = input.Mean;
	double readDistance = 0;

	if (l.Position + lLen >= lRefLen || r.Position + rLen >= rRefLen)
		return;
	if (leftAlg.MapQuality < input.MapQ || rightAlg.MapQuality < input.MapQ)
		return;
	if (lLen < input.MinReadLength || rLen < input.MinReadLength)
		return;
	int leftEdit, rightEdit;
	if (!leftAlg.GetTag("NM", leftEdit) || !rightAlg.GetTag("NM", rightEdit) || leftEdit > input.MaxEditDistance || rightEdit > input.MaxEditDistance)
		return;

	if (!input.IsIllumina)
	{
		if (!l.IsReverseStrand && !r.IsReverseStrand)
		{
			distance += lRefLen - l.Position - lLen + r.Position;
			readDistance += rRefLen - r.Position + l.Position + lLen;
		}
		else if (!l.IsReverseStrand && r.IsReverseStrand)
		{
			distance += lRefLen - l.Position - lLen + rRefLen - r.Position - rLen;
			readDistance += r.Position + rLen + l.Position + lLen;
		}
		else if (l.IsReverseStrand && !r.IsReverseStrand)
		{
			distance += l.Position + r.Position;
			readDistance += rRefLen - r.Position + lRefLen - l.Position;
		}
		else if (l.IsReverseStrand && r.IsReverseStrand)
		{
			distance += l.Position + rRefLen - r.Position - rLen;
			readDistance += r.Position + rLen + lRefLen - l.Position;
		}
	}
	else
	{
		if (!l.IsReverseStrand && !r.IsReverseStrand)
		{
			distance += l.Position + r.Position;
			readDistance += lRefLen - l.Position + rRefLen - r.Position;
		}
		else if (!l.IsReverseStrand && r.IsReverseStrand)
		{
			distance += l.Position + rRefLen - r.Position - rLen;
			readDistance += lRefLen - l.Position + r.Position + rLen;
		}
		else if (l.IsReverseStrand && !r.IsReverseStrand)
		{
			distance += lRefLen - l.Position - lLen + r.Position;
			readDistance += lLen + l.Position + rRefLen - r.Position;
		}
		else if (l.IsReverseStrand && r.IsReverseStrand)
		{
			distance += lRefLen - l.Position - lLen + rRefLen - r.Position - rLen;
			readDistance += lLen + l.Position + r.Position + rLen;
		}
	}

	if (noOverlapDeviation > Helpers::Eps && readDistance > input.Mean + noOverlapDeviation * input.Std)
		return;

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

