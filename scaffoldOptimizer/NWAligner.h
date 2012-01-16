/*
 * scaffoldOptimizer : solves the MIQP optimization and produces linear scaffold
 * sequences.
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

#ifndef _NWALIGNER_H
#define	_NWALIGNER_H

#include <memory>
#include <algo/align/nw/mm_aligner.hpp>
#include "Sequence.h"

using namespace std;

class NWAligner
{
public:
    NWAligner(const FastASequence &a, const FastASequence &b, int matchScore = 2, int mismatchScore = -3, const string &endSpace = ((char *)"xxxx"));
    
public:
    int Align();
    int GetAlignmentScore();
    FastASequence GetAlignmentA() const;
    FastASequence GetAlignmentB() const;
    FastASequence GetConsensus() const;
    const FastASequence &GetSequenceA() const;
    const FastASequence &GetSequenceB() const;

private:
    string getFormattedAlignment() const;
    
private:
    bool aligned;
    int score;
    auto_ptr<ncbi::CNWAligner> aligner;
    FastASequence seqA, seqB;
};

#endif	/* _NWALIGNER_H */
