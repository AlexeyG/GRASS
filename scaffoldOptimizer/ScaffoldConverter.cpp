#include "ScaffoldConverter.h"
#include <string>
#include "Helpers.h"
#include "NWAligner.h"
//#include "MinMax.h"
#include "OverlapperConfiguration.h"
#include "ContigOverlapper.h"

#include <iostream>

using namespace std;

vector<FastASequence> ScaffoldConverter::ToFasta(const DataStore &store, const Scaffold &scaffold, const OverlapperConfiguration &config)
{
    vector<FastASequence> ans;
    string name, sequence;
    int count = scaffold.ContigCount();
    int solutionEnd = 0;
    int actualEnd = 0;
    int scaffoldOffset = 0;

    for (int i = 0; i < count; i++)
    {
        ScaffoldContig contig = scaffold[i];
        FastASequence contigSeq = store[contig.Id].Sequence;
        string sign;
        int contigLen = contigSeq.Nucleotides.length();
        int contigEnd = (!contig.T ? contig.X + contigLen : contig.X);
        if (i == 0)
        {
            scaffoldOffset = (!contig.T ? contig.X : contig.X - contigLen + 1);
        }
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
        double overlapScore;
        int overlapOffset = ContigOverlapper::FindBestOverlap(sequence, contigSeq.Nucleotides, solutionDistance, config, overlapScore, overlapConsensus);
        cout << "Got offset: " << overlapOffset << " of Q: " << overlapScore << endl;
        if (overlapScore * 100 > config.OverlapQuality) // we have a good overlap
        {
            //
        }
        else // we don't have a good overlap
        {
            if (solutionDistance >= 0) // and we were not supposed to overlap
            {
                string spacer(solutionDistance, 'N');
                sequence = sequence + spacer + contigSeq.Nucleotides;
                actualEnd += solutionDistance + contigLen;
            }
            else if (solutionDistance < config.NoSplitOverlapLength) // the predicted overlap is short
            {
                sequence = sequence + contigSeq.Nucleotides;
                actualEnd += contigLen;
            }
            else // we predicted a long overlap
            {
                // split
                ans.push_back(FastASequence(sequence, name));
                sequence = contigSeq.Nucleotides;
                name.clear();
                // update coords
                actualEnd = contigLen;
                scaffoldOffset = (!contig.T ? contig.X : contig.X - contigLen + 1);
            }
            name = (!name.empty() ? name + "|" + sign : sign) + Helpers::ItoStr(contig.Id);
            solutionEnd = max(solutionEnd, contigEnd - scaffoldOffset);
        }
    }

    if (!sequence.empty())
        ans.push_back(FastASequence(sequence, name));
    
    return ans;
}

vector<FastASequence> ScaffoldConverter::ToFasta(const DataStore &store, const vector<Scaffold> &scaffold, const OverlapperConfiguration &config)
{
    int count = scaffold.size();
    vector<FastASequence> seq;
    for (int i = 0; i < count; i++)
    {
        vector<FastASequence> tmp = ToFasta(store, scaffold[i], config);
        seq.insert(seq.end(), tmp.begin(), tmp.end());
    }

    return seq;
}
