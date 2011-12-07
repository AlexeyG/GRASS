#include "BreakpointCount.h"

#include <iostream>

using namespace std;

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
        cout << "Orientation: " << a.QueryID << endl;
        cout << "             " << a.ReferencePosition << " - " << b.ReferencePosition << endl;
        return true;
    }
    if ((a.ReferencePosition >= b.ReferencePosition) && (a.IsReferenceReverse == a.IsQueryReverse))
    {
        Total++, Order++;
        cout << "Order: " << a.QueryID << endl;
        cout << "       " << a.ReferencePosition << " - " << b.ReferencePosition << endl;
        return true;
    }
    int distanceReference = b.ReferencePosition - a.ReferencePosition + a.ReferenceAlignmentLength;
    int distanceQuery = b.QueryPosition - a.QueryPosition + a.QueryAlignmentLength;
    if (isDistanceBreakpoint(distanceReference, distanceQuery))
    {
        Total++, Distance++;
        cout << "Distance: " << a.QueryID << endl;
        cout << "       " << a.ReferencePosition << " - " << b.ReferencePosition << endl;
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
        i = j;
        j++;
    }
    return count;
}

bool BreakpointCount::isDistanceBreakpoint(int distA, int distB)
{
    return abs(distA - distB) > DistanceThreshold;
}