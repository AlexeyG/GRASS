#include "MummerAligner.h"
#include "Defines.h"
#include <cstdio>

MummerAligner::MummerAligner()
{
    alignmentString = "";
}

bool MummerAligner::Align(const string& referenceFileName, const string& scaffoldsFileName)
{
    FILE *pipe = popen(cmd, "r");
    if (!pipe)
        return false;
    char buffer[MAX_LINE];
    alignmentString = "";
    while(!feof(pipe))
    {
        if(fgets(buffer, MAX_LINE, pipe) != NULL)
            alignmentString += buffer;
        else
            break;
    }
    return !alignmentString.empty();
}
