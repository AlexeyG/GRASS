#include "ScaffoldConverter.h"
#include <string>
#include "Helpers.h"
#include "NWAligner.h"
#include "MinMax.h"

#include <iostream>

using namespace std;

FastASequence ScaffoldConverter::ToFasta(const DataStore &store, const Scaffold &scaffold)
{
    string name, sequence;
    int count = scaffold.ContigCount();
    int solutionEnd = 0;
    int actualEnd = 0;
    for (int i = 0; i < count; i++)
    {
        ScaffoldContig contig = scaffold[i];
        FastASequence contigSeq = store[contig.Id].Sequence;
        int contigLen = contigSeq.Nucleotides.length();
        if (!contig.T)
        {
            int solutionDistance = contig.X - solutionEnd;
            string spacer(max(solutionDistance, 0), 'N'); // spacer sequence
            if (solutionDistance < 0)
            {
                int overlapLength = -solutionDistance;
                FastASequence seqA = FastASequence(sequence.substr(actualEnd - overlapLength, overlapLength), "prev");
                FastASequence seqB = FastASequence(contigSeq.Nucleotides.substr(0, overlapLength), "current");
                NWAligner aligner(seqA, seqB);
                int score = aligner.Align();
                cout << "Overlap: " << overlapLength << " | " << score << endl;
            }
            sequence = sequence + spacer + contigSeq.Nucleotides;
            actualEnd += spacer.length() + contigLen;
            solutionEnd += spacer.length() + contigLen;
        }
        else
        {
            contigSeq.ReverseCompelement();
            
            int solutionDistance = contig.X - contigLen + 1 - solutionEnd;
            string spacer(max(solutionDistance, 0), 'N'); // spacer sequence
            if (solutionDistance < 0)
            {
                int overlapLength = -solutionDistance;
                FastASequence seqA = FastASequence(sequence.substr(actualEnd - overlapLength, overlapLength), "prev");
                FastASequence seqB = FastASequence(contigSeq.Nucleotides.substr(0, overlapLength), "current");
                NWAligner aligner(seqA, seqB);
                int score = aligner.Align();
                cout << "Overlap: " << overlapLength << " | " << score << endl;
            }
            
            sequence = sequence + spacer + contigSeq.Nucleotides;
            actualEnd += spacer.length() + contigLen;
            solutionEnd += spacer.length() + contigLen;
        }
        string sign = (!contig.T ? "-" : "+");
        name = (!name.empty() ? name + "|" + sign : sign) + Helpers::ItoStr(contig.Id);
    }

    return FastASequence(sequence, name);
}

vector<FastASequence> ScaffoldConverter::ToFasta(const DataStore &store, const vector<Scaffold> &scaffold)
{
    int count = scaffold.size();
    vector<FastASequence> seq(count);
    for (int i = 0; i < count; i++)
        seq[i] = ToFasta(store, scaffold[i]);

    return seq;
}
