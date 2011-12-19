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
        
#ifdef _OVERLAP
        
        string overlapConsensus;
        double overlapScore;
        int overlapOffset = ContigOverlapper::FindBestOverlap(sequence, contigSeq.Nucleotides, solutionDistance, config, overlapScore, overlapConsensus);
        cout << "Got offset: " << overlapOffset << " of Q: " << overlapScore << " | Predicted: " << solutionDistance << endl;
        int leftLen = actualEnd, rightLen = contigLen, overlapLen = 0;
        if (overlapScore * 100 > config.OverlapQuality) // we have a good overlap
        {
            if (leftLen <= rightLen)
            {
                /*
                 *        ----  L
                 * ----------   R
                 */
                if (overlapOffset < -rightLen) // right overhang
                {
                    overlapLen = leftLen + rightLen + overlapOffset; // 4 + 10 - 1199
                    sequence = contigSeq.Nucleotides.substr(0, contigLen - overlapLen) + overlapConsensus + sequence.substr(overlapLen, actualEnd - overlapLen);
                }
                /*
                 * ----        L
                 * ----------  R
                 */
                else if (overlapOffset <= -leftLen) // left inside
                {
                    overlapLen = leftLen;
                    sequence = contigSeq.Nucleotides.substr(0, -overlapOffset - leftLen) + overlapConsensus + contigSeq.Nucleotides.substr(-overlapOffset, contigLen + overlapOffset);
                }
                /*
                 * ----         L
                 *  ----------  R
                 */
                else // left overhang
                {
                    overlapLen = -overlapOffset;
                    sequence = sequence.substr(0, actualEnd - overlapLen) + overlapConsensus + contigSeq.Nucleotides.substr(overlapLen, contigLen - overlapLen);
                }
            }
            /*
             * ----------
             *         -----
             */
            else
            {
                /*
                 * ----------   L
                 *        ----  R
                 */
                if (overlapOffset > -rightLen) // right overhang
                {
                    overlapLen = -overlapOffset;
                    sequence = sequence.substr(0, actualEnd - overlapLen) + overlapConsensus + contigSeq.Nucleotides.substr(overlapLen, contigLen - overlapLen);
                }
                /*
                 * ----------  L
                 * ----        R
                 */
                else if (overlapOffset > -leftLen) // right inside
                {
                    overlapLen = rightLen;;
                    // leftSequence = left.substr(leftLen + offset, overlapLen);
                    sequence = sequence.substr(0, actualEnd + overlapOffset) + overlapConsensus + sequence.substr(-overlapOffset, -overlapOffset - contigLen);
                }
                /*
                 *  ----------  L
                 * ----         R
                 */
                else // left overhang
                {
                    overlapLen = leftLen + rightLen + overlapOffset;
                    sequence = contigSeq.Nucleotides.substr(0, contigLen - overlapLen) + overlapConsensus + sequence.substr(overlapLen, actualEnd - overlapLen);
                }
            }
        }
        else // we don't have a good overlap
        {
            if (solutionDistance >= 0) // and we were not supposed to overlap
            {
                string spacer(solutionDistance, 'N');
                sequence = sequence + spacer + contigSeq.Nucleotides;
            }
            else if (-solutionDistance <= config.NoSplitOverlapLength) // the predicted overlap is short
            {
                sequence = sequence + contigSeq.Nucleotides;
            }
            else // we predicted a long overlap
            {
                cout << "SPLITTING!" << endl;
                // split
                ans.push_back(FastASequence(sequence, name));
                sequence = contigSeq.Nucleotides;
                name.clear();
                // update coords
                scaffoldOffset = (!contig.T ? contig.X : contig.X - contigLen + 1);
            }
        }
        name = (!name.empty() ? name + "|" + sign : sign) + Helpers::ItoStr(contig.Id);
        actualEnd = sequence.length();
        solutionEnd = max(solutionEnd, contigEnd - scaffoldOffset);
    }

#else
    for (int i = 0; i < count; i++)
    {
        ScaffoldContig contig = scaffold[i];
        FastASequence contigSeq = store[contig.Id].Sequence;
        int contigLen = contigSeq.Nucleotides.length();
        if (!contig.T)
        {
            //cout << i << " Forward: " << contig.X << " - " << end << " = " << contig.X - end << endl;
            string spacer(max(contig.X - solutionEnd, 0), 'N');
            sequence = sequence + spacer + contigSeq.Nucleotides;
            name = (!name.empty() ? name + "|+" : "+") + Helpers::ItoStr(contig.Id);
            solutionEnd += spacer.length() + contigLen;
        }
        else
        {
            //cout << i << " Reverse: " << contig.X << " - " << contigLen << " + 1 -" << end << " = " << contig.X - contigLen + 1 - end << endl;
            contigSeq.ReverseCompelement();
            string spacer(max(contig.X - contigLen + 1 - solutionEnd, 0), 'N');
            sequence = sequence + spacer + contigSeq.Nucleotides;
            name = (!name.empty() ? name + "|-" : "-") + Helpers::ItoStr(contig.Id);
            solutionEnd += spacer.length() + contigLen;
        }
    }
    
#endif
    
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
