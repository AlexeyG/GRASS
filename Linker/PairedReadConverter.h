#ifndef _PAIREDREADCONVERTER_H
#define _PAIREDREADCONVERTER_H

#include "Configuration.h"
#include "DataStore.h"
#include "XATag.h"

class PairedReadConverter
{
public:
	PairedReadConverter(DataStore &store);
	static bool IsCorrectRelativeOrientation(const XATag &l, const XATag &r, bool isIllumina);
	enum PairedReadConverterResult { Success, FailedLeftAlignment, FailedRightAlignment, FailedLeftConversion, FailedRightConversion, FailedLinkCreation };

public:
	PairedReadConverterResult Process(const Configuration &config, const PairedInput &input);

private:
	PairedReadConverterResult alignAndConvert(const Configuration &config, const PairedInput &input);
	PairedReadConverterResult createLinksFromAlignment(int groupId, int maxHits, const PairedInput &input);
	void createLinksForPair(int groupId, const BamAlignment &leftAlg, const vector<XATag> &leftTags, const BamAlignment &rightAlg, const vector<XATag> &rightTags, const PairedInput &input, int maxHits);
	void addLinkForTagPair(int groupId, const XATag &l, int lLen, const XATag &r, int rLen, const PairedInput &input, int factor = 1);
	void removeBamFiles();

private:
	DataStore &dataStore;
	string leftBamFileName;
	string rightBamFileName;
};
#endif
