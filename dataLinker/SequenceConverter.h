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
    enum SequenceConverterResult { Success, FailedAlignment, FailedReadAlignment, FailedReadSequences, FailedLinkCreation, InconsistentReferenceSets };

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
