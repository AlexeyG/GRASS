#include "XATag.h"

XATag::XATag(const BamAlignment &alg)
{
	RefID = alg.RefID;
	Position = alg.Position;
	IsReverseStrand = alg.IsReverseStrand();
}
