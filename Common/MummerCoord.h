/* 
 * File:   MummerCoord.h
 * Author: alexeyg
 *
 * Created on December 6, 2011, 9:35 AM
 */

#ifndef _MUMMERCOORD_H
#define	_MUMMERCOORD_H

class MummerCoord
{
public:
    MummerCoord(int referenceID = 0, int queryID = 0, int referencePosition = 0, int queryPosition = 0, int referenceAlignmentLength = 0, int queryAlignmentLength = 0, bool isReferenceReverse = false, bool isQueryReverse = false, double identity = 1)
    : ReferenceID(referenceID), QueryID(queryID), ReferencePosition(referencePosition), QueryPosition(queryPosition),
      ReferenceAlignmentLength(referenceAlignmentLength), QueryAlignmentLength(queryAlignmentLength),
      IsReferenceReverse(isReferenceReverse), IsQueryReverse(isQueryReverse) {};
    
public:
    bool operator< (const MummerCoord &b) const
    {
        if (this->QueryID < b.QueryID)
            return true;
        if (this->QueryID > b.QueryID)
            return false;
        if (this->QueryPosition < b.QueryPosition)
            return true;
        return false;
    }
      
public:
    // 0-based reference sequence ID
    int ReferenceID;
    // 0-based query sequence ID
    int QueryID;
    // 0-based position in reference sequence
    int ReferencePosition;
    // 0-based position in query sequence
    int QueryPosition;
    // Reference alignment length
    int ReferenceAlignmentLength;
    // Query alignment length
    int QueryAlignmentLength;
    bool IsReferenceReverse;
    bool IsQueryReverse;
    // [0; 1] identity percentage for this alignment
    double Identity;
};

#endif	/* _MUMMERCOORD_H */
