/* 
 * File:   ReadCoverage.h
 * Author: alexeyg
 *
 * Created on November 29, 2011, 3:41 PM
 */

#ifndef _READCOVERAGE_H
#define	_READCOVERAGE_H

#include <vector>

using namespace std;

class ReadCoverage
{
public:
    ReadCoverage(int contigCount = 0);
    
public:
    int GetContigCount() const;
    void SetContigCount(int count);
    void AddLocation(int id, int location);
    void UpdateAverage(int readLength);
    void SetAverageReadLength(long long totalReadLength, int totalReadCount);
    
public:
    vector< vector<int> > ReadLocations;
    double AverageReadLength;
    int TotalReadCount;
    long long TotalReadLength;
};

#endif	/* READCOVERAGE_H */
