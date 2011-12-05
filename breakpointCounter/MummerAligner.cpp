#include "MummerAligner.h"
#include "Defines.h"
#include <cstdio>

MummerAligner::MummerAligner()
{
    alignmentString = "";
}

bool MummerAligner::Align(const string& referenceFileName, const string& scaffoldsFileName)
{
    char buffer[MAX_LINE];
    alignmentString = "";
    sprintf(buffer, "nucmer ", referenceFileName.c_str(), scaffoldsFileName.c_str());
    FILE *pipe = popen(buffer, "r");
    if (pipe == NULL)
        return false;
    while(!feof(pipe))
    {
        if(fgets(buffer, MAX_LINE, pipe) != NULL)
            alignmentString += buffer;
        else
            break;
    }
    return !alignmentString.empty();
}
