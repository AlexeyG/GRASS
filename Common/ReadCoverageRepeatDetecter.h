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


/* 
 * File:   ReadCoverageRepeatDetecter.h
 * Author: alexeyg
 *
 * Created on November 30, 2011, 9:28 AM
 */


#ifndef _READCOVERAGEREPEATDETECTER_H
#define	_READCOVERAGEREPEATDETECTER_H

#include <vector>
#include "ReadCoverage.h"
#include "DataStore.h"

using namespace std;

class ReadCoverageRepeatDetecter
{   
public:
    static vector<int> Detect(double expectedCoverage, const ReadCoverage &coverage, const DataStore &store, double uniquenessCutoff);
};

#endif	/* _READCOVERAGEREPEATDETECTER_H */
