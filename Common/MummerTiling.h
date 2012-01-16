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


#ifndef _MUMMERTILING_H
#define	_MUMMERTILING_H

class MummerTiling
{
public:
    MummerTiling(int referenceID = 0, int queryID = 0, int referencePosition = 0, int referenceLength = 0, int queryLength = 0, bool isReverse = false, double identity = 1, double coverage = 1)
    : ReferenceID(referenceID), QueryID(queryID), ReferencePosition(referencePosition),
      ReferenceLength(referenceLength), QueryLength(queryLength),
      IsReverse(isReverse), Identity(identity), Coverage(coverage) {};
    
public:
    bool operator< (const MummerTiling &b) const
    {
        if (ReferenceID < b.ReferenceID)
            return true;
        if (ReferenceID > b.ReferenceID)
            return false;
        if (ReferencePosition < b.ReferencePosition)
            return true;
        if (ReferencePosition > b.ReferencePosition)
            return false;
        if (ReferenceLength > b.ReferenceLength)
            return true;
        if (ReferenceLength < b.ReferenceLength)
            return false;
        if (QueryLength > b.QueryLength)
            return true;
        if (QueryLength < b.QueryLength)
            return false;
        if (QueryID < b.QueryID)
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
    // Reference alignment length
    int ReferenceLength;
    // Query length
    int QueryLength;
    // Is the query reverse?
    bool IsReverse;
    // [0; 1] identity percentage for this alignment
    double Identity;
    // [0; 1] coverage percentage for this alignment
    double Coverage;
};

#endif

