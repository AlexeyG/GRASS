#ifndef _PAIREDREADPROCESSOR_H
#define _PAIREDREADPROCESSOR_H

#include "Configuration.h"
#include "DataStore.h"
#include "XATag.h"

class PairedReadProcessor
{
public:
	PairedReadProcessor();
	static bool IsCorrectRelativeOrientation(const XATag &l, const XATag &r, bool isIllumina);
	enum PairedReadProcessorResult { Success, FailedLeftAlignment, FailedRightAlignment, FailedLeftConversion, FailedRightConversion, FailedIO };

public:
	PairedReadProcessorResult Process(const Configuration &config, const PairedInput &input);

private:
	PairedReadProcessorResult alignAndConvert(const Configuration &config, const PairedInput &input);
	PairedReadProcessorResult processAlignment(const PairedInput &input);
	void removeBamFiles();

private:
	string leftBamFileName;
	string rightBamFileName;
};
#endif
