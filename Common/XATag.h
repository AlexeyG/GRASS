/*
 * Common : a collection of classes (re)used throughout the scaffolder implementation.
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


#ifndef _XATAG_H
#define _XATAG_H

#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include <string>

using namespace std;
using namespace BamTools;

class XATag
{
public:
	XATag(int position, int refID, bool isReverseStrand) : Position(position), RefID(refID), IsReverseStrand(isReverseStrand) {};
	XATag(const BamAlignment &alg);
	XATag(const string &str, const BamReader &reader);

public:
	int Position;
	int RefID;
	bool IsReverseStrand;
};
#endif
