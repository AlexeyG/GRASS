#include "ContigOverlapper.h"
#include "NWAligner.h"

#include <iostream>

using namespace std;

int ContigOverlapper::FindEndOverlap(const string &left, const string &right, int distance, const OverlapperConfiguration &config, int &score, string &consensus)
{
    if (left.empty() || right.empty()) // empty sequences do not overlap
        return 0;
    
    if (distance == -13464)
    {
        cout << "Traced" << endl;
    }
    
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
//        cout << "Left: " << left << " -> " << leftSequence << endl;
//        cout << "Right: " << right << " -> " << rightSequence << endl;
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
        
        /*cout << "Got score: " << score << endl;
        cout << "Got length: " << offset << " : " << consensus << endl;
        cout << leftAlignment << endl;
        cout << rightAlignment << endl;*/
        
        return offset;
    }
}
