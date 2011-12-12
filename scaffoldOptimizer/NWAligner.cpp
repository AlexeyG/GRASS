#include "NWAligner.h"
#include <algo/align/nw/nw_formatter.hpp>
#include <string>
#include <sstream>

using namespace std;

NWAligner::NWAligner(const FastASequence &a, const FastASequence &b, int matchScore, int mismatchScore, const string &endSpace)
    : seqA(a), seqB(b)
{
    aligned = false;
    score = 0;
    aligner = auto_ptr<ncbi::CNWAligner>(new ncbi::CMMAligner(a.Nucleotides, b.Nucleotides));
    aligner->SetWm(matchScore);
    aligner->SetWms(mismatchScore);
    aligner->SetScoreMatrix(NULL);
    if (endSpace.length() == 4)
        aligner->SetEndSpaceFree(endSpace[0] == 'z', endSpace[1] == 'z', endSpace[2] == 'z', endSpace[3] == 'z');
}

int NWAligner::Align()
{
    if (aligned)
        return score;
    aligned = true;
    score = aligner->Run();
    cout << "Got score: " << score << endl;
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
    getline(ss, alignment); // skip header 2
    getline(ss, alignment);
    
    return FastASequence(alignment, "NW|" + seqB.Name());
}

FastASequence NWAligner::GetConsensus() const
{
    if (!aligned)
        return FastASequence();
    
    FastASequence a = GetSequenceA(), b = GetSequenceB();
    int length = a.Nucleotides.length();
    string str("N", length);
    for (int i = 0; i < length; i++)
    {
        if (a.Nucleotides[i] == b.Nucleotides[i])
            str[i] = a.Nucleotides[i];
        else if (a.Nucleotides[i] == '-' || a.Nucleotides[i] == 'N' || a.Nucleotides[i] == 'n')
            str[i] = b.Nucleotides[i];
        else if (b.Nucleotides[i] == '-' || b.Nucleotides[i] == 'N' || b.Nucleotides[i] == 'n')
            str[i] = a.Nucleotides[i];
        else
            str[i] = 'N'; // could to a better job here
    }
    
    return FastASequence(str, a.Name() + "|" + b.Name());
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
