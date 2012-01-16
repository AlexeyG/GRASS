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

#include "PairedReadProcessor.h"
#include "Aligner.h"
#include "Converter.h"
#include "Helpers.h"
#include "AlignmentReader.h"
#include "Writer.h"
#include <sstream>

PairedReadProcessor::PairedReadProcessor()
{
}

bool PairedReadProcessor::IsCorrectRelativeOrientation(const XATag &l, const XATag &r, bool isIllumina)
{
	return (isIllumina ? l.IsReverseStrand != r.IsReverseStrand : l.IsReverseStrand == r.IsReverseStrand);
}

PairedReadProcessor::PairedReadProcessorResult PairedReadProcessor::Process(const Configuration &config, const PairedInput &input)
{
	PairedReadProcessorResult result = alignAndConvert(config, input);
	if (result == Success)
		result = processAlignment(input);
	removeBamFiles();
	return result;
}

PairedReadProcessor::PairedReadProcessorResult PairedReadProcessor::alignAndConvert(const Configuration &config, const PairedInput &input)
{
	Aligner *leftAlignment, *rightAlignment;
	PairedReadProcessorResult result = Success;
	bool success = true;

	leftAlignment = (input.IsIllumina ? (Aligner *)new BWAAligner(config.ReferenceFileName, input.LeftFileName, config.BWAConfig) : new NovoAlignAligner(config.ReferenceFileName, input.LeftFileName, config.NovoAlignConfig));
	rightAlignment = (input.IsIllumina ? (Aligner *)new BWAAligner(config.ReferenceFileName, input.RightFileName, config.BWAConfig) : new NovoAlignAligner(config.ReferenceFileName, input.RightFileName, config.NovoAlignConfig));

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

PairedReadProcessor::PairedReadProcessorResult PairedReadProcessor::processAlignment(const PairedInput &input)
{
	PairedReadProcessorResult result = Success;
	AlignmentReader leftReader, rightReader;
	FastQWriter leftWriter, rightWriter;
	if (!leftReader.Open(leftBamFileName) || !rightReader.Open(rightBamFileName))
		result = FailedIO;
	if (!leftWriter.Open(input.OutputPrefix + "_1.fastq") || !rightWriter.Open(input.OutputPrefix + "_2.fastq"))
		result = FailedIO;
	vector<XATag> leftTags, rightTags;
	BamAlignment leftAlignment, rightAlignment;
	while (leftReader.GetNextAlignmentGroup(leftAlignment, leftTags) && rightReader.GetNextAlignmentGroup(rightAlignment, rightTags))
	{
		int leftCount = leftTags.size(), rightCount = rightTags.size();
		/*if (leftAlignment.Name.find("96138.1") != string::npos || leftAlignment.Name.find("96138.2") != string::npos)
		{
			cout << "left: " << leftAlignment.Name << " : " << leftCount << endl;
			cout << "right: " << rightAlignment.Name << " : " << rightCount << endl;
		}*/
		if (leftCount > 1 || rightCount > 1)
			continue;
		if (!leftWriter.Write(FastQSequence(leftAlignment)) || !rightWriter.Write(FastQSequence(rightAlignment)))
		{
			result = FailedIO;
			break;
		}
	}
	leftReader.Close();
	rightReader.Close();
	leftWriter.Close();
	rightWriter.Close();
	return result;
}

void PairedReadProcessor::removeBamFiles()
{
	if (!leftBamFileName.empty())
		Helpers::RemoveFile(leftBamFileName);
	if (!rightBamFileName.empty())
		Helpers::RemoveFile(rightBamFileName);
	leftBamFileName.clear();
	rightBamFileName.clear();
}
