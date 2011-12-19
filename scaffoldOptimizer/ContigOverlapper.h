/* 
 * File:   ContigOverlapper.h
 * Author: alexeyg
 *
 * Created on December 11, 2011, 3:39 PM
 */

#ifndef _CONTIGOVERLAPPER_H
#define	_CONTIGOVERLAPPER_H

#include "OverlapperConfiguration.h"
#include <string>

using namespace std;

class ContigOverlapper
{
public:
    static int FindBestOverlap(const std::string &left, const std::string &right, int distance, const OverlapperConfiguration &config, double &score, std::string &consensus);
    static int GetOverlapSequences(int offset, const std::string &left, const std::string &right, std::string &leftSequence, std::string &rightSequence);
    
private:
    static int GetAlignmentScore(const std::string &left, const std::string &right);
    static std::string GetConsensus(const std::string &left, const std::string &right);
};

#endif	/* _CONTIGOVERLAPPER_H */
