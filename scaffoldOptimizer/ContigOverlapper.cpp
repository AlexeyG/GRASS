#include "ContigOverlapper.h"
//#include "MinMax.h"
#include "NWAligner.h"

int ContigOverlapper::FindEndOverlap(const string &left, const string &right, int distance, const OverlapperConfiguration &config, int &score, string &consensus)
{
    if (distance >= 0 && distance - config.InitialOverlapDeviation >= 0) // Sequences don't overlap even with offset
        return 0;
    else
    {
        int lenLeft = left.length();
        int lenRight = right.length();
        int offset = min(min(lenLeft, lenRight), config.InitialOverlapDeviation - distance);
        if (offset - config.InitialOverlapDeviation > config.MaximumAlignmentLength) // Too long to perform global alignment
            return -1;
        
        string leftSequence = left.substr(left.length() - offset, offset);
        string rightSequence = right.substr(0, offset);
        NWAligner aligner(FastASequence(leftSequence, "left"), FastASequence(rightSequence, "right"));
        aligner.Align();
        
        string leftAlignment = aligner.GetAlignmentA().Nucleotides;
        string rightAlignment = aligner.GetAlignmentB().Nucleotides;
        int leftPointer = leftAlignment.length() - 1;
        int rightPointer = 0;
        while (leftPointer >= rightPointer && leftAlignment[leftPointer] == '-' && rightAlignment[rightPointer] == '-')
        {
            leftPointer--;
            rightPointer++;
        }
        offset = leftPointer - rightPointer + 1;
        if (offset <= 0) // we should not actually overlap
            return 0;
        
        leftSequence = left.substr(left.length() - offset, offset);
        rightSequence = right.substr(0, offset);
        NWAligner aligner2(FastASequence(leftSequence, "left"), FastASequence(rightSequence, "right"));
        score = aligner2.Align();
        consensus = aligner2.GetConsensus().Nucleotides;
        return offset;
    }
}
