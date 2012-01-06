#include "SequenceConverter.h"
#include "Aligner.h"
#include "MummerCoordReader.h"
#include "Reader.h"
#include <algorithm>

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
    cout << "Got " << coords.size() << " alignments!" << endl;
    for (int i = 0; i < (int)coords.size(); i++)
    {
        cout << coords[i].ReferenceID << "\t" << coords[i].ReferencePosition << "\t" << coords[i].QueryID << "\t" << coords[i].ReferenceAlignmentLength << endl;
    }
    
    return result;
}

int SequenceConverter::addLinkGroup(const SequenceInput &input)
{
    LinkGroup group("Reference sequence alignment", "Alignment to sequences " + input.FileName " with " << input.Std << " of weight " << input.Weight);
    return dataStore.AddGroup(group);
}

SequenceConverter::SequenceConverterResult SequenceConverter::alignContigs(const string &sequenceFileName, const Configuration &config, Coords &coords)
{
    MummerAligner aligner(sequenceFileName, config.InputFileName, config.MummerConfig);
    std::vector<FastASequence> references;
    FastAReader faReader;
    if (!faReader.Open(sequenceFileName) || faReader.Read(references) <= 0)
        return FailedReadSequences;
    MummerCoordReader reader(references, dataStore.GetContigs());
    if (!aligner.Align())
        return FailedAlignment;
    if (!reader.Open(aligner.OutputFileName) || reader.Read(coords) == 0)
        return FailedReadAlignment;
    return Success;
}

void SequenceConverter::createLinksFromAlignment(int groupID, Coords &coords, const SequenceInput &input)
{
    auto prev = coords.begin();
    for (auto cur = prev + 1; cur != coords.end(); prev++, cur++)
    {
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
        }
    }
}

int SequenceConverter::filterAlignmentsOnLength(Coords &coords, double minBases)
{
    int count = 0;
    Coords newCoords;
    for (auto it = coords.begin(); it != coords.end(); it++)
        if (it->ReferenceAlignmentLength > minBases && it->QueryAlignmentLength > minBases)
            newCoords.push_back(*it);
        else
            count++;
    coords.swap(newCoords);
    return count;
}

void SequenceConverter::sortAlignments(Coords &coords)
{
    sort(coords.begin(), coords.end(), referencePositionComparator);
}

bool SequenceConverter::referencePositionComparator(const MummerCoord &l, const MummerCoord &r)
{
    if (l.ReferenceID < r.ReferenceID)
        return true;
    if (l.ReferenceID > r.ReferenceID)
        return false;
    if (l.ReferencePosition < r.ReferencePosition)
        return true;
    if (l.ReferencePosition > r.ReferencePosition)
        return false;
    if (l.ReferenceAlignmentLength > r.ReferenceAlignmentLength)
        return true;
    if (l.ReferenceAlignmentLength < r.ReferenceAlignmentLength)
        return false;
    if (l.QueryAlignmentLength > r.QueryAlignmentLength)
        return true;
    if (l.QueryAlignmentLength < r.QueryAlignmentLength)
        return false;
    if (l.QueryID < r.QueryID)
        return true;
    return false;
}
