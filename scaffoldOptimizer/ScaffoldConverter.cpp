#include "ScaffoldConverter.h"
#include <string>
#include "Helpers.h"
#include "NWAligner.h"
//#include "MinMax.h"
#include "OverlapperConfiguration.h"
#include "ContigOverlapper.h"

#include <iostream>

using namespace std;

FastASequence ScaffoldConverter::ToFasta(const DataStore &store, const Scaffold &scaffold, const OverlapperConfiguration &config)
{
    string name, sequence;
    int count = scaffold.ContigCount();
    int solutionEnd = 0;
    int actualEnd = 0;
    for (int i = 0; i < count; i++)
    {
        ScaffoldContig contig = scaffold[i];
        FastASequence contigSeq = store[contig.Id].Sequence;
        string sign;
        int contigLen = contigSeq.Nucleotides.length();
        int solutionDistance;
        if (!contig.T)
        {
            solutionDistance = contig.X - solutionEnd;
            sign = "+";
        }
        else
        {
            contigSeq.ReverseCompelement();
            solutionDistance = contig.X - contigLen + 1 - solutionEnd;
            sign = "-";
        }
        
        string overlapConsensus;
        int overlapScore;
        bool shouldBreak = false;
        int overlapLength = ContigOverlapper::FindEndOverlap(sequence, contigSeq.Nucleotides, solutionDistance, config, overlapScore, overlapConsensus);
        if (overlapLength == 0) // should not overlap
        {
            if (solutionDistance >= 0) // and we were not supposed to overlap
            {
                string spacer(solutionDistance, 'N'); // spacer sequence
                sequence = sequence + spacer + contigSeq.Nucleotides;
                actualEnd += spacer.length() + contigLen;
                solutionEnd += spacer.length() + contigLen;
            }
            else if (-solutionDistance <= config.NoSplitOverlapLength) // we were supposed to overlap, but the overlap was not too long
            {
                sequence = sequence + contigSeq.Nucleotides;
                actualEnd += contigLen;
                solutionEnd += solutionDistance + contigLen;
            }
            else if (!sequence.empty()) // predicted overlap was too long!
            {
                cout << "Predicted overlap of " << -solutionDistance << " but didn't find it -> have to split!" << endl;
                cout << "Seq: " << sequence << endl;
                shouldBreak = true;
            }
        }
        else if (overlapLength > 0) // should overlap
        {
            double overlapQuality = (double)overlapScore / (double)(overlapLength * 2.0) * 100;
            if (overlapQuality >= config.OverlapQuality) // found overlap is of high quality
            {
                sequence = sequence.substr(0, sequence.length() - overlapLength) + overlapConsensus + contigSeq.Nucleotides.substr(overlapLength, contigLen - overlapLength);
                actualEnd += contigLen - 2 * overlapLength + overlapConsensus.length();
                solutionEnd += contigLen + solutionDistance;
            }
            else if (solutionDistance >= 0) // overlap is of poor quality, but we didn't want to overlap anyways
            {
                string spacer(solutionDistance, 'N');
                sequence = sequence + spacer + contigSeq.Nucleotides;
                actualEnd += spacer.length() + contigLen;
                solutionEnd += spacer.length() + contigLen;
            }
            else if (-solutionDistance <= config.NoSplitOverlapLength) // supposed to overlap, but the predicted overlap was not too long
            {
                sequence = sequence + contigSeq.Nucleotides;
                actualEnd += contigLen;
                solutionEnd += solutionDistance + contigLen;
            }
            else
            {
                cout << "Predicted overlap of " << -solutionDistance << " but found overlap of " << overlapLength << "bp of quality " << overlapQuality << " -> have to split!" << endl;
                shouldBreak = true;
            }
        }
        else // did not perform alignment
        {
            if (solutionDistance >= 0 || -solutionDistance <= config.NoSplitOverlapLength)
            {
                string spacer(max(solutionDistance, 0), 'N'); // spacer sequence
                sequence = sequence + spacer + contigSeq.Nucleotides;
                actualEnd += spacer.length() + contigLen;
                solutionEnd += contigLen + solutionDistance;
            }
            else
            {
                cout << "Predicted overlap of " << -solutionDistance << " but could not perform alignment -> have to split!" << endl;
                shouldBreak = true;
                // should probably split? Ask what Dick thinks about it.
            }
        }
        
        if (!shouldBreak)
            name = (!name.empty() ? name + "|" + sign : sign) + Helpers::ItoStr(contig.Id);
        else
        {
            // handle scaffold breaking here!
        }
    }

    return FastASequence(sequence, name);
}

vector<FastASequence> ScaffoldConverter::ToFasta(const DataStore &store, const vector<Scaffold> &scaffold, const OverlapperConfiguration &config)
{
    int count = scaffold.size();
    vector<FastASequence> seq(count);
    for (int i = 0; i < count; i++)
        seq[i] = ToFasta(store, scaffold[i], config);

    return seq;
}
