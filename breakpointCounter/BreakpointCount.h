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

using namespace std;

class BreakpointCount
{
public:
    BreakpointCount(int distanceThreshold = 1000) : DistanceThreshold(distanceThreshold), Joins(0), Order(0), Orientation(0), Distance(0), Total(0) {}
    
public:
    bool IsBreakpoint(const MummerCoord &a, const MummerCoord &b);
    int ProcessAlignments(const vector<MummerCoord> &coords, const vector<FastASequence> &references);
    
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
    
public:
    static void Sort(vector<MummerCoord> &coords);
};

#endif	/* _BREAKPOINTCOUNT_H */
