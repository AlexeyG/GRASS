#include "Sequence.h"
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;

string FastASequence::Name() const
{
	string copy(Comment);
	stringstream in(copy);
	string name;
	in >> name;
	return name;
}

Sequence::Sequence(const BamAlignment &alg)
{
	Nucleotides = alg.QueryBases;
	if (alg.IsReverseStrand())
	{
		reverse(Nucleotides.begin(), Nucleotides.end());
		complement();
	}
}

void Sequence::ReverseCompelement()
{
	reverse(Nucleotides.begin(), Nucleotides.end());
	complement();
}

void Sequence::complement()
{
	int n = Nucleotides.length();
	for (int i = 0; i < n; i++)
		switch (Nucleotides[i])
		{
			case 'A' : Nucleotides[i] = 'T';
				break;
			case 'T' : Nucleotides[i] = 'A';
				break;
			case 'C' : Nucleotides[i] = 'G';
				break;
			case 'G' : Nucleotides[i] = 'C';
				break;
			case 'a' : Nucleotides[i] = 't';
				break;
			case 't' : Nucleotides[i] = 'a';
				break;
			case 'c' : Nucleotides[i] = 'g';
				break;
			case 'g' : Nucleotides[i] = 'c';
				break;
			case 'N' : Nucleotides[i] = 'N';
				break;
			case 'n' : Nucleotides[i] = 'n';
				break;
			default :
				Nucleotides[i] = 'N';
		}
}

FastASequence::FastASequence(const BamAlignment &alg)
	: Sequence(alg), Comment(alg.Name)
{
}

FastQSequence::FastQSequence(const BamAlignment &alg)
	: FastASequence(alg), Quality(alg.Qualities)
{
	if (alg.IsReverseStrand())
		reverse(Quality.begin(), Quality.end());
}

void FastQSequence::ReverseCompelement()
{
	reverse(Nucleotides.begin(), Nucleotides.end());
	complement();
	reverse(Quality.begin(), Quality.end());
}
