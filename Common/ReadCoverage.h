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
    ReadCoverage(int contigCount);
    
public:
    int GetContigCount() const;
    void SetContigCount(int count);
    void AddLocation(int id, int location);
    void UpdateAverage(int readLength);
    
public:
    vector< vector<int> > ReadLocations;
    double AverageReadLength;
    int TotalReadCount;
    
private:
    long long totalReadLength;
};

#endif	/* READCOVERAGE_H */
