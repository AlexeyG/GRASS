#include "ContigOverlapper.h"
#include "NWAligner.h"
#include "Configuration.h"

#include <iostream>
#include <cstdio>

using namespace std;

int ContigOverlapper::FindBestOverlap(const std::string &left, const std::string &right, int predictedGap, const OverlapperConfiguration &config, double &score, std::string &consensus)
{
    if (left.empty() || right.empty()) // empty sequences do not overlap
        return 0;
    
    int leftLen = left.length(), rightLen = right.length();
    
    int maxDeviation = config.InitialOverlapDeviation;
    int maxOffset = predictedGap + maxDeviation;
    int minOffset = predictedGap - maxDeviation;
    
    if (minOffset >= 0) // the furthest we are allowed to go does not overlap
    {
        score = -1;
        return predictedGap;
    }
    
    if (maxOffset >= 0)
        maxOffset = -1;
    if (minOffset <= -(leftLen + rightLen))
        minOffset = -(leftLen + rightLen) + 1;
    
    double bestScore = -1;
    int bestOffset = 0;
    
    for (int offset = minOffset; offset <= maxOffset; offset++)
    {
        double distanceScore = (maxDeviation - abs(predictedGap - offset)) / (double) maxDeviation;
        std::string leftSequence, rightSequence;
        int overlapLen = GetOverlapSequences(offset, left, right, leftSequence, rightSequence);
        // check if overlapLen is too large to align!
        if (overlapLen <= config.MaximumAlignmentLength)
        {
            double alignmentScore = (double)GetAlignmentScore(leftSequence, rightSequence) / (double)(2 * overlapLen); // [-1.5; 1] alignment score
            double score = alignmentScore * distanceScore;
            if (score > bestScore)
            {
                bestScore = score;
                bestOffset = offset;
            }
            //printf("Offset: %i Score: %.6lf (%i)\n", offset, score, overlapLen);
        }
    }
    
    if (bestScore != -1) // actually should get consensus
    {
        std::string leftSequence, rightSequence;
        GetOverlapSequences(bestOffset, left, right, leftSequence, rightSequence);
        consensus = GetConsensus(leftSequence, rightSequence);
    }
    score = bestScore;
    return bestOffset;
}

int ContigOverlapper::GetOverlapSequences(int offset, const std::string &left, const std::string &right, std::string &leftSequence, std::string &rightSequence)
{
    int leftLen = left.length(), rightLen = right.length();
    leftSequence.clear(); rightSequence.clear();
    int overlapLen = 0;
    
    if (offset >= 0)
        return 0;
    
    /*
     * ----
     *   ----------
     */
    if (leftLen <= rightLen)
    {
        if (offset <= -(leftLen + rightLen))
        {
            overlapLen = 0;
        }
        /*
         *        ----  L
         * ----------   R
         */
        else if (offset < -rightLen) // right overhang
        {
            overlapLen = leftLen + rightLen + offset; // 4 + 10 - 1199
            leftSequence = left.substr(0, overlapLen);
            rightSequence = right.substr(rightLen - overlapLen, overlapLen);
        }
        /*
         * ----        L
         * ----------  R
         */
        else if (offset <= -leftLen) // left inside
        {
            overlapLen = leftLen;
            leftSequence = left;
            rightSequence = right.substr(-offset - leftLen, leftLen);
        }
        /*
         * ----         L
         *  ----------  R
         */
        else // left overhang
        {
            overlapLen = -offset;
            leftSequence = left.substr(leftLen + offset, overlapLen);
            rightSequence = right.substr(0, overlapLen);
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
        if (offset > -rightLen) // right overhang
        {
            overlapLen = -offset;
            leftSequence = left.substr(leftLen + offset, overlapLen);
            rightSequence = right.substr(0, overlapLen);
        }
        /*
         * ----------  L
         * ----        R
         */
        else if (offset > -leftLen) // right inside
        {
            overlapLen = rightLen;
            leftSequence = left.substr(leftLen + offset, overlapLen);
            rightSequence = right;
        }
        /*
         *  ----------  L
         * ----         R
         */
        else if (offset > -(leftLen + rightLen)) // left overhang
        {
            overlapLen = leftLen + rightLen + offset;
            leftSequence = left.substr(0, overlapLen);
            rightSequence = right.substr(rightLen - overlapLen, overlapLen);
        }
        else
        {
            overlapLen = 0;
        }
    }
    return overlapLen;
}

int ContigOverlapper::GetAlignmentScore(const std::string &left, const std::string &right)
{
    //cout << left << endl;
    //cout << right << endl;
    //cout << "-----" << endl;
    NWAligner aligner(FastASequence(left, "l"), FastASequence(right, "r"));
    return aligner.Align();
}

std::string ContigOverlapper::GetConsensus(const std::string &left, const std::string &right)
{
    NWAligner aligner(FastASequence(left, "l"), FastASequence(right, "r"));
    aligner.Align();
    return aligner.GetConsensus().Nucleotides;
}
