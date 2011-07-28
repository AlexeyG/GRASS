#include "XATag.h"

XATag::XATag(const BamAlignment &alg)
{
	RefID = alg.RefID;
	Position = alg.Position;
	IsReverseStrand = alg.IsReverseStrand();
}

XATag::XATag(const string &str, const BamReader &reader)
{
        if (str.empty())
                throw new exception();
        int n = str.length();
        int i = 0, s = 0;
        while (i < n - 1 && str[i] != ',') i++;
        string name = str.substr(0, i);
        s = ++i; while (i < n - 1 && str[i] != ',') i++;
        string position = str.substr(s, i - s);
        if (position.empty())
                throw new exception();
        if (position[0] != '+' && position[0] != '-')
                throw new exception();
        RefID = reader.GetReferenceID(name);
        Position = atoi(position.c_str() + 1);
        IsReverseStrand = position[0] == '-';
}
