#include "BreakpointCount.h"
#include "Configuration.h"

//#include <iostream>
#include <algorithm>

using namespace std;

BreakpointCount::BreakpointCount(int distanceThreshold)
        : DistanceThreshold(distanceThreshold), Joins(0), Order(0), Orientation(0), Distance(0), Total(0)
{
    scaffoldCoverage = auto_ptr<Depth>(new Depth());
    referenceCoverage = auto_ptr<Depth>(new Depth());
}

BreakpointCount::~BreakpointCount()
{
}

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
        //cout << "Orientation: " << a.QueryID << endl;
        //cout << "             " << a.ReferencePosition << " - " << b.ReferencePosition << endl;
        return true;
    }
    if ((a.ReferencePosition >= b.ReferencePosition) && (a.IsReferenceReverse == a.IsQueryReverse))
    {
        Total++, Order++;
        //cout << "Order: " << a.QueryID << endl;
        //cout << "       " << a.ReferencePosition << " - " << b.ReferencePosition << endl;
        return true;
    }
    int distanceReference = getQueryDistance(a, b);
    int distanceQuery = getReferenceDistance(a,b);
    if (isDistanceBreakpoint(distanceReference, distanceQuery))
    {
        Total++, Distance++;
        //cout << "Distance: " << a.QueryID << endl;
        //cout << "(" << abs(distanceReference - distanceQuery) << "): " << distanceReference << " - " << distanceQuery << endl;
        //cout << "       " << a.ReferencePosition << " - " << b.ReferencePosition << endl;
        return true;
    }
    return false;
}

int BreakpointCount::ProcessAlignments(const vector<MummerCoord> &coords, const vector<FastASequence> &references, const vector<FastASequence> &scaffolds)
{
    resizeCoverageVector(*referenceCoverage, references);
    resizeCoverageVector(*scaffoldCoverage, scaffolds);
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

double BreakpointCount::GetReferenceCoverage() const
{
    long long total = 0;
    long long covered = 0;
    int nSeq = referenceCoverage->size();
    for (int i = 0; i < nSeq; i++)
    {
        int seqLen = referenceCoverage->at(i).size();
        total += seqLen;
        for (int j = 0; j < seqLen; j++)
            if ((referenceCoverage->at(i))[j])
                covered++;
    }
    return (double)covered / (double)total;
}

double BreakpointCount::GetScaffoldCoverage() const
{
    long long total = 0;
    long long covered = 0;
    int nSeq = scaffoldCoverage->size();
    for (int i = 0; i < nSeq; i++)
    {
        int seqLen = scaffoldCoverage->at(i).size();
        total += seqLen;
        for (int j = 0; j < seqLen; j++)
            if ((scaffoldCoverage->at(i))[j])
                covered++;
    }
    return (double)covered / (double)total;
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

bool BreakpointCount::isDistanceBreakpoint(int distA, int distB) const
{
    return abs(distA - distB) > DistanceThreshold;
}

void BreakpointCount::coverSequences(const MummerCoord &c)
{
    if (c.ReferenceID >= (int)referenceCoverage->size())
        return;
    if (c.QueryID >= (int)scaffoldCoverage->size())
        return;
    
    for (int i = 0; i < c.ReferenceAlignmentLength; i++)
        (referenceCoverage->at(c.ReferenceID))[c.ReferencePosition + i] = true;
    for (int i = 0; i < c.QueryAlignmentLength; i++)
        (scaffoldCoverage->at(c.QueryID)[c.QueryPosition + i] = true;
}

void BreakpointCount::Sort(vector<MummerCoord> &coords)
{
    sort(coords.begin(), coords.end());
}

int BreakpointCount::getQueryDistance(const MummerCoord &a, const MummerCoord &b)
{
    return b.QueryPosition - a.QueryPosition - a.QueryAlignmentLength;
}

int BreakpointCount::getReferenceDistance(const MummerCoord &a, const MummerCoord &b)
{
    return (a.IsQueryReverse == a.IsReferenceReverse ? b.ReferencePosition - a.ReferencePosition - a.ReferenceAlignmentLength : a.ReferencePosition - b.ReferencePosition - b.ReferenceAlignmentLength);
}

void BreakpointCount::resizeCoverageVector(Depth &depth, const vector<FastASequence> &seq)
{
    int nSeq = seq.size();
    depth.resize(nSeq, vector<bool>());
    for (int i = 0; i < nSeq; i++)
        depth[i].resize(seq[i].Nucleotides.length());
}
