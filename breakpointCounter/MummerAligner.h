/* 
 * File:   MummerAligner.h
 * Author: alexeyg
 *
 * Created on December 5, 2011, 9:31 PM
 */

#ifndef _MUMMERALIGNER_H
#define	_MUMMERALIGNER_H

#include <string>

using namespace std;

class MummerAligner
{
public:
    MummerAligner();
    
public:
    bool Align(const string &referenceFileName, const string &scaffoldsFileName);
    bool ParseAlignment();
    
private:
    string alignmentString;
};

#endif	/* _MUMMERALIGNER_H */
