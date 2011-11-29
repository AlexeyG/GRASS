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
    ReadCoverage ReadCoverage;
        
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
