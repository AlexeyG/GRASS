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

#include "ScaffoldConverter.h"
#include <string>
#include "Helpers.h"
#include "NWAligner.h"
//#include "MinMax.h"
#include "OverlapperConfiguration.h"
#include "ContigOverlapper.h"

#include <iostream>

//#define _OVERLAP

using namespace std;

vector<FastASequence> ScaffoldConverter::ToFasta(const DataStore &store, const Scaffold &scaffold, const OverlapperConfiguration &config)
{
    vector<FastASequence> ans;
    string name, sequence;
    int count = scaffold.ContigCount();
    int solutionEnd = 0;

#ifdef _OVERLAP
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
        //cout << "Got offset: " << overlapOffset << " of Q: " << overlapScore << " | Predicted: " << solutionDistance << endl;
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
                    //cout << "A: " << leftLen << " + " << rightLen << " = " << sequence.length() << endl;
                }
                /*
                 * ----        L
                 * ----------  R
                 */
                else if (overlapOffset <= -leftLen) // left inside
                {
                    overlapLen = leftLen;
                    sequence = contigSeq.Nucleotides.substr(0, -overlapOffset - leftLen) + overlapConsensus + contigSeq.Nucleotides.substr(-overlapOffset, contigLen + overlapOffset);
                    //cout << "B: " << leftLen << " + " << rightLen << " = " << sequence.length() << endl;
                }
                /*
                 * ----         L
                 *  ----------  R
                 */
                else // left overhang
                {
                    overlapLen = -overlapOffset;
                    sequence = sequence.substr(0, actualEnd - overlapLen) + overlapConsensus + contigSeq.Nucleotides.substr(overlapLen, contigLen - overlapLen);
                    //cout << "C: " << leftLen << " + " << rightLen << " = " << sequence.length() << endl;
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
                    //cout << "D: " << leftLen << " + " << rightLen << " = " << sequence.length() << endl;
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
                    //cout << "E: " << leftLen << " + " << rightLen << " = " << sequence.length() << endl;
                }
                /*
                 *  ----------  L
                 * ----         R
                 */
                else // left overhang
                {
                    overlapLen = leftLen + rightLen + overlapOffset;
                    sequence = contigSeq.Nucleotides.substr(0, contigLen - overlapLen) + overlapConsensus + sequence.substr(overlapLen, actualEnd - overlapLen);
                    //cout << "F: " << leftLen << " + " << rightLen << " = " << sequence.length() << endl;
                }
            }
        }
        else // we don't have a good overlap
        {
            if (solutionDistance >= 0) // and we were not supposed to overlap
            {
                string spacer(solutionDistance, 'N');
                sequence = sequence + spacer + contigSeq.Nucleotides;
                //cout << "G: " << leftLen << " + " << rightLen << " = " << sequence.length() << " (" << solutionDistance << ")" << endl;
            }
            else if (-solutionDistance <= config.NoSplitOverlapLength) // the predicted overlap is short
            {
                sequence = sequence + contigSeq.Nucleotides;
                //cout << "H: " << leftLen << " + " << rightLen << " = " << sequence.length() << endl;
            }
            else // we predicted a long overlap
            {
                //cout << "SPLITTING!" << endl;
                // split
                ans.push_back(FastASequence(sequence, name));
                sequence = contigSeq.Nucleotides;
                name.clear();
                // update coords
                scaffoldOffset = (!contig.T ? contig.X : contig.X - contigLen + 1);
                solutionEnd = 0;
                //cout << "I: " << leftLen << " + " << rightLen << " = " << sequence.length() << endl;
            }
        }
        //cout << "Left: " << sequence << endl;
        name = (!name.empty() ? name + "|" + sign : sign) + Helpers::ItoStr(contig.Id);
        actualEnd = sequence.length();
        //cout << "max of " << solutionEnd << " and " << contigEnd << endl;
        solutionEnd = max(solutionEnd, contigEnd);
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
            string spacer(max((int)contig.X - solutionEnd, 0), 'N');
            sequence = sequence + spacer + contigSeq.Nucleotides;
            name = (!name.empty() ? name + "|+" : "+") + Helpers::ItoStr(contig.Id);
            solutionEnd += spacer.length() + contigLen;
        }
        else
        {
            //cout << i << " Reverse: " << contig.X << " - " << contigLen << " + 1 -" << end << " = " << contig.X - contigLen + 1 - end << endl;
            contigSeq.ReverseCompelement();
            string spacer(max((int)contig.X - contigLen + 1 - solutionEnd, 0), 'N');
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
