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
    static int FindEndOverlap(const std::string &left, const std::string &right, int distance, const OverlapperConfiguration &config, int &score, std::string &consensus);
};

#endif	/* _CONTIGOVERLAPPER_H */
