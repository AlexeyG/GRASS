#include "NWAligner.h"
#include <algo/align/nw/nw_formatter.hpp>
#include <string>

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
    aligner->GetAlignment();
    return score;
}

int NWAligner::GetAlignmentScore()
{
    return Align();
}

FastASequence NWAligner::GetAlignment() const
{
    if (!aligned)
        return FastASequence();
    string alignment;
    ncbi::CNWFormatter formatter(*aligner);
    formatter.AsText(&alignment, ncbi::CNWFormatter::eFormatFastA, seqA.Nucleotides.length() + seqB.Nucleotides.length());
    cout << "Alg: " << alignment << endl;
    return FastASequence(alignment, "NW|" + seqA.Name() + "|" + seqB.Name());
}

const FastASequence &NWAligner::GetSequenceA() const
{
    return seqA;
}

const FastASequence &NWAligner::GetSequenceB() const
{
    return seqB;
}
