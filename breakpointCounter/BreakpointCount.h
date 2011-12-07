/* 
 * File:   BreakpointCount.h
 * Author: alexeyg
 *
 * Created on December 6, 2011, 8:51 PM
 */

#ifndef _BREAKPOINTCOUNT_H
#define	_BREAKPOINTCOUNT_H

#include "MummerCoord.h"
#include "Sequence.h"
#include <vector>
#include <memory>

using namespace std;

typedef vector< vector<bool> > Depth;

class BreakpointCount
{
public:
    BreakpointCount(int distanceThreshold = 10000);
    ~BreakpointCount();
    
public:
    bool IsBreakpoint(const MummerCoord &a, const MummerCoord &b);
    int ProcessAlignments(const vector<MummerCoord> &coords, const vector<FastASequence> &references, const vector<FastASequence> &scaffolds);
    double GetReferenceCoverage() const;
    double GetScaffoldsCoverage() const;
    
public:
    static void Sort(vector<MummerCoord> &coords);
    
public:
    int DistanceThreshold;
    
public:
    int Joins;
    int Order;
    int Orientation;
    int Distance;
    int Total;
    
private:
    int processAlignmentGroup(vector<MummerCoord>::const_iterator start, vector<MummerCoord>::const_iterator finish);
    bool isDistanceBreakpoint(int distA, int distB) const;
    void coverSequences(const MummerCoord &c);
    
private:
    static int getQueryDistance(const MummerCoord &a, const MummerCoord &b);
    static int getReferenceDistance(const MummerCoord &a, const MummerCoord &b);
    static void resizeCoverageVector(Depth &depth, const vector<FastASequence> &seq);
    
private:
    auto_ptr<Depth> scaffoldCoverage;
    auto_ptr<Depth> referenceCoverage;
};

#endif	/* _BREAKPOINTCOUNT_H */
