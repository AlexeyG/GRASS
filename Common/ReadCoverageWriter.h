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
 * File:   ReadCoverageWriter.h
 * Author: alexeyg
 *
 * Created on November 29, 2011, 4:17 PM
 */

#ifndef _READCOVERAGEWRITER_H
#define	_READCOVERAGEWRITER_H

#include <cstddef>
#include <string>
#include "ReadCoverage.h"

using namespace std;

class ReadCoverageWriter {
public:
    ReadCoverageWriter();
    virtual ~ReadCoverageWriter();

public:
    bool Open(const string &fileName, const string &mode = "wb");
    bool Close();
    bool Write(const ReadCoverage &coverage);

protected:
	FILE *out;
};

#endif	/* _READCOVERAGEWRITER_H */

