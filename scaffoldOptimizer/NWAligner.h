/* 
 * File:   NWAligner.h
 * Author: alexeyg
 *
 * Created on December 9, 2011, 2:18 AM
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
    NWAligner(const FastASequence &a, const FastASequence &b, int matchScore = 2, int mismatchScore = -3, const string &endSpace = ((char *)"xzzx"));
    
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
