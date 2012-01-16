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
