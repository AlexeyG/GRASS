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

#ifndef _SEQUENCECONVERTER_H
#define	_SEQUENCECONVERTER_H

#include "Configuration.h"
#include "DataStore.h"
#include "MummerCoord.h"

#include <string>
#include <vector>

typedef std::vector<MummerCoord> Coords;

class SequenceConverter
{
public:
    SequenceConverter(DataStore &store);
    enum SequenceConverterResult { Success, FailedAlignment, FailedReadAlignment, FailedReadSequences };

public:
    SequenceConverterResult Process(const Configuration &config, const SequenceInput &input);

private:
    SequenceConverterResult alignContigs(const std::string &sequenceFileName, const Configuration &config, Coords &coords);
    int addLinkGroup(const SequenceInput &input);
    void createLinksFromAlignment(int groupID, Coords &coords, const SequenceInput &input);
    static int filterAlignmentsOnLength(Coords &coords, double minBases);
    static void sortAlignments(Coords &coords);
    static bool referencePositionComparator(const MummerCoord &l, const MummerCoord &r);
    
private:
    DataStore &dataStore;
};

#endif	/* _SEQUENCECONVERTER_H */
