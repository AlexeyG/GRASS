/*
 * dataLinker : creates abstract contig links from the available information sources.
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

#ifndef _PAIREDREADCONVERTER_H
#define _PAIREDREADCONVERTER_H

#include "Configuration.h"
#include "DataStore.h"
#include "XATag.h"
#include "ReadCoverage.h"
#include <vector>

using namespace std;

class PairedReadConverter
{
public:
	PairedReadConverter(DataStore &store);
	static bool IsCorrectRelativeOrientation(const XATag &l, const XATag &r, bool isIllumina);
        enum PairedReadConverterResult { Success, FailedLeftAlignment, FailedRightAlignment, FailedLeftConversion, FailedRightConversion, FailedLinkCreation, InconsistentReferenceSets };

public:
	PairedReadConverterResult Process(const Configuration &config, const PairedInput &input);

public:
    ReadCoverage ContigReadCoverage;
        
private:
	PairedReadConverterResult alignAndConvert(const Configuration &config, const PairedInput &input);
	PairedReadConverterResult createLinksFromAlignment(int groupId, int maxHits, const PairedInput &input, double noOverlapDeviation);
        void createLinksForPair(int groupId, const BamAlignment &leftAlg, const vector<XATag> &leftTags, const BamAlignment &rightAlg, const vector<XATag> &rightTags, const PairedInput &input, double noOverlapDeviation, int maxHits);
        void processCoverage(const BamAlignment &alg, const vector<XATag> &tags);
	void addLinkForTagPair(int groupId, const XATag &l, const BamAlignment &leftAlg, const XATag &r, const BamAlignment &rightAlg, const PairedInput &input, double noOverlapDeviation, int factor = 1);
	void removeBamFiles();

private:
	DataStore &dataStore;
	string leftBamFileName;
	string rightBamFileName;
};
#endif
