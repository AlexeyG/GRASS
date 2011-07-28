#include "ContigInfo.h"

ContigInfo::ContigInfo(const FastASequence &contig)
	: Contig(contig)
{
	Length = Contig.Nucleotides.length();
	Position = determineContigPosition();
}

bool ContigInfo::operator< (const ContigInfo &other) const
{
	return Position < other.Position;
}

int ContigInfo::determineContigPosition()
{
	if (Length == 0)
		return -1;
	int i = 0;
	if (Contig.Comment[i] == '+' || Contig.Comment[i] == '-')
		i++;
	int j = i + 1;
	while (j < Length && Contig.Comment[j] != '|')
		j++;
	if (j >= Length)
		return -1;
	int num = 0;
	for (; i < j; i++)
		num = num * 10 + Contig.Comment[i] - '0';
	if (Contig.Comment[0] == '-')
		num -= Length;
	return num;
}
