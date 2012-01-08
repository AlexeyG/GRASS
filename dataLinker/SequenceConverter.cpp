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

#include "SequenceConverter.h"
#include "Aligner.h"
#include "MummerTilingReader.h"
#include "Reader.h"
#include <algorithm>
#include <sstream>

using namespace std;

SequenceConverter::SequenceConverter(DataStore &store)
     : dataStore(store) 
{
}

SequenceConverter::SequenceConverterResult SequenceConverter::Process(const Configuration &config, const SequenceInput &input)
{
    SequenceConverterResult result;
    Coords coords;
    result = alignContigs(input.FileName, config, coords);
    if (result != Success)
        return result;
    filterAlignmentsOnLength(coords, input.MinAlignmentLength);
    sortAlignments(coords);
    createLinksFromAlignment(addLinkGroup(input), coords, input);
    
    return result;
}

int SequenceConverter::addLinkGroup(const SequenceInput &input)
{
    stringstream ss;
    ss << "Alignment to sequences " << input.FileName << " with " << input.Std << " of weight " << input.Weight;
    LinkGroup group("Reference sequence alignment", ss.str());
    return dataStore.AddGroup(group);
}

SequenceConverter::SequenceConverterResult SequenceConverter::alignContigs(const string &sequenceFileName, const Configuration &config, Coords &coords)
{
    MummerTiler tiler(sequenceFileName, config.InputFileName, config.MummerConfig);
    std::vector<FastASequence> references;
    FastAReader faReader;
    if (!faReader.Open(sequenceFileName) || faReader.Read(references) <= 0)
        return FailedReadSequences;
    MummerTilingReader reader(references, dataStore.GetContigs());
    if (!tiler.Align())
        return FailedAlignment;
    if (!reader.Open(tiler.OutputFileName) || reader.Read(coords) == 0)
        return FailedReadAlignment;
    return Success;
}

void SequenceConverter::createLinksFromAlignment(int groupID, Coords &coords, const SequenceInput &input)
{
    cout << "Got " << coords.size() << " alignments!" << endl;
    auto prev = coords.begin();
    for (auto cur = prev + 1; cur != coords.end(); prev++, cur++)
    {
        cout << prev->ReferenceID << "\t" << prev->ReferencePosition << "\t" << prev->QueryID << "\t" << prev->ReferenceAlignmentLength << endl;
        if (prev->ReferenceID == cur->ReferenceID && prev->QueryID != cur->QueryID)
        {
            int prevLength = dataStore[prev->QueryID].Sequence.Nucleotides.length();
            int curLength = dataStore[cur->QueryID].Sequence.Nucleotides.length();
            int distance = prev->QueryPosition + (cur->ReferencePosition - prev->ReferencePosition + 1) + (curLength - cur->QueryPosition - cur->QueryAlignmentLength + 1);
            bool equalOrientation = (prev->IsQueryReverse ^ prev->IsReferenceReverse) == (cur->IsQueryReverse ^ cur->IsReferenceReverse);
            bool forwardOrder = prev->IsQueryReverse == prev->IsReferenceReverse;
            double weight = input.Weight * cur->Identity * cur->QueryAlignmentLength / (double) curLength * prev->QueryAlignmentLength / (double) prevLength;
            ContigLink link(prev->QueryID, cur->QueryID, distance, input.Std, equalOrientation, forwardOrder, weight);
            dataStore.AddLink(groupID, link);
            cout << "   added link between " << prev->QueryID << " and " << cur->QueryID << endl;
        }
    }
}

void SequenceConverter::sortAlignments(Coords &coords)
{
    sort(coords.begin(), coords.end());
}
