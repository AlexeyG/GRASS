/* 
 * File:   OverlapperConfiguration.h
 * Author: alexeyg
 *
 * Created on December 11, 2011, 1:34 PM
 */

#ifndef _OVERLAPPERCONFIGURATION_H
#define	_OVERLAPPERCONFIGURATION_H

class OverlapperConfiguration
{
public:
    OverlapperConfiguration();
    
public:
    int MaximumAlignmentLength;
    int InitialOverlapDeviation;
    int NoSplitOverlapLength;
    double OverlapQuality;
};

#endif	/* _OVERLAPPERCONFIGURATION_H */
