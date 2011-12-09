#include "NWAligner.h"
#include <algo/align/nw/nw_formatter.hpp>
#include <string>

using namespace std;

NWAligner::NWAligner(const FastASequence &a, const FastASequence &b)
    : seqA(a), seqB(b)
{
    aligned = false;
    score = 0;
    aligner = auto_ptr<ncbi::CNWAligner>(ncbi::CMMAligner(a.Nucleotides, b.Nucleotides, &NCBISM_Blosum62));
}

int NWAligner::Align() const
{
    if (aligned)
        return score;
    aligned = true;
    return score = aligner->Run();
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
    formatter.AsText(&alignment, CNWFormatter::eFormatFastA, seqA.Nucleotides.length() + seqB.Nucleotides.length());
    return FastASequence(alignment, "Algn|" + seqA.Name() + "|" + seqB.Name());
}

const FastASequence &NWAligner::GetSequenceA() const
{
    return seqA;
}

const FastASequence &NWAligner::GetSequenceB() const
{
    return seqB;
}
