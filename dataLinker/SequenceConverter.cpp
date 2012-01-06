#include "SequenceConverter.h"
#include "Aligner.h"
#include "MummerCoordReader.h"
#include "Reader.h"

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
    cout << "Got " << coords.size() << " alignments!" << endl;
    
    return result;
}

SequenceConverter::SequenceConverterResult SequenceConverter::alignContigs(const string &sequenceFileName, const Configuration &config, Coords &coords)
{
    MummerAligner aligner(sequenceFileName, config.InputFileName, config.MummerConfig);
    std::vector<FastASequence> references;
    FastAReader faReader;
    if (!faReader.Open(config.InputFileName) || faReader.Read(references) <= 0)
        return FailedReadSequences;
    MummerCoordReader reader(references, dataStore.GetContigs());
    if (!aligner.Align())
        return FailedAlignment;
    if (!reader.Open(aligner.OutputFileName) || reader.Read(coords) == 0)
        return FailedReadAlignment;
    return Success;
}

int SequenceConverter::filterAlignmentsOnLength(Coords &coords, double minBases)
{
    int count = 0;
    Coords newCoords;
    for (auto it = coords.begin(); it != coords.end(); it++)
        if (it->ReferenceAlignmentLength > minBases && it->QueryAlignmentLenght > minBases)
            newCoords.push_back(*it);
        else
            count++;
    coords.swap(newCoords);
    return count;
}

/*int SequenceConverter::filterAlignments(Coords &coords, double minBases)
{
    int count = 0;
    Coords newCoords;
    for (auto it = coords.begin(); it != coords.end(); it++)
        if (it->Identity * it->ReferenceAlignmentLength > minBases)
            newCoords.push_back(*it);
        else
            count++;
    coords.swap(newCoords);
    return count;
}*/
