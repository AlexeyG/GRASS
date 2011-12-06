#include "BreakpointCount.h"

bool BreakpointCount::IsBreakpoint(const MummerCoord &a, const MummerCoord &b)
{
    if (a.ReferenceID != b.ReferenceID)
    {
        Total++, Joins++;
        return true;
    }
    if ((a.IsReferenceReverse ^ b.IsReferenceReverse) != (a.IsQueryReverse ^ b.IsQueryReverse))
    {
        Total++, Orientation++;
        return true;
    }
    if ((a.ReferencePosition >= b.ReferencePosition) && (a.IsReferenceReverse == a.IsQueryReverse))
    {
        Total++, Order++;
        return true;
    }
    return false;
}

int BreakpointCount::ProcessAlignments(const vector<MummerCoord> &coords, const vector<FastASequence> &references)
{
    int count = 0;
    vector<MummerCoord>::const_iterator it = coords.begin();
    while (it != coords.end())
    {
        auto start = it;
        while (it != coords.end() && it->QueryID == start->QueryID)
            it++;
        count += processAlignmentGroup(start, it);
    }
    return count;
}

int BreakpointCount::processAlignmentGroup(vector<MummerCoord>::const_iterator start, vector<MummerCoord>::const_iterator finish)
{
    int count = 0;
    auto i = start;
    auto j = start + 1;
    while (i != finish)
    {
        // mark coverage
        if (j != finish)
            if (IsBreakpoint(*i, *j))
                count++;
    }
    return count;
}
