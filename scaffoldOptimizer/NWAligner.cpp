#include "NWAligner.h"
#include <algo/align/nw/nw_formatter.hpp>
#include <string>
#include <sstream>

using namespace std;

NWAligner::NWAligner(const FastASequence &a, const FastASequence &b, int matchScore, int mismatchScore)
    : seqA(a), seqB(b)
{
    aligned = false;
    score = 0;
    //&NCBISM_Blosum62
    aligner = auto_ptr<ncbi::CNWAligner>(new ncbi::CMMAligner(a.Nucleotides, b.Nucleotides));
    aligner->SetWm(matchScore);
    aligner->SetWms(mismatchScore);
    aligner->SetScoreMatrix(NULL);
}

int NWAligner::Align()
{
    if (aligned)
        return score;
    aligned = true;
    score = aligner->Run();
    cout << GetAlignmentA().Nucleotides << endl;
    cout << GetAlignmentB().Nucleotides << endl;
    return score;
}

int NWAligner::GetAlignmentScore()
{
    return Align();
}

FastASequence NWAligner::GetAlignmentA() const
{
    if (!aligned)
        return FastASequence();
    string alignment;
    stringstream ss(getFormattedAlignment());
    getline(ss, alignment); // skip header 1
    getline(ss, alignment);
    
//    cout << "Alg: " << alignment << endl;
    return FastASequence(alignment, "NW|" + seqA.Name());
}

FastASequence NWAligner::GetAlignmentB() const
{
    if (!aligned)
        return FastASequence();
    string alignment;
    stringstream ss(getFormattedAlignment());
    getline(ss, alignment); // skip header 1
    getline(ss, alignment); // skip sequence 1
    getline(ss, alignment); // skup header 2
    getline(ss, alignment);
    
//    cout << "Alg: " << alignment << endl;
    return FastASequence(alignment, "NW|" + seqB.Name());
}

FastASequence NWAligner::GetConsensusAlignment() const
{
    if (!aligned)
        return FastASequence();
    // calculate consensus here
    return FastASequence();
}

const FastASequence &NWAligner::GetSequenceA() const
{
    return seqA;
}

const FastASequence &NWAligner::GetSequenceB() const
{
    return seqB;
}

string NWAligner::getFormattedAlignment() const
{
    string alignment;
    ncbi::CNWFormatter formatter(*aligner);
    formatter.AsText(&alignment, ncbi::CNWFormatter::eFormatFastA, seqA.Nucleotides.length() + seqB.Nucleotides.length());
    return alignment;
}
