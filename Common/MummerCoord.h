/*
 * Common : a collection of classes (re)used throughout the scaffolder implementation.
 * Copyright (C) 2011  Alexey Gritsenko
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 * 
 * 
 * Email: a.gritsenko@tudelft.nl
 * Mail: Delft University of Technology
 *       Faculty of Electrical Engineering, Mathematics, and Computer Science
 *       Department of Mediamatics
 *       P.O. Box 5031
 *       2600 GA, Delft, The Netherlands
 */


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
      IsReferenceReverse(isReferenceReverse), IsQueryReverse(isQueryReverse), Identity(identity) {};
    
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
