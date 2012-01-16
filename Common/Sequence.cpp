/*
 * Common : a collection of classes (re)used throughout the scaffolder implementation.
 * Copyright (C) 2011  Alexey Gritsenko
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 * 
 * 
 * Email: a.gritsenko@tudelft.nl
 * Mail: Delft University of Technology
 *       Faculty of Electrical Engineering, Mathematics, and Computer Science
 *       Department of Mediamatics
 *       P.O. Box 5031
 *       2600 GA, Delft, The Netherlands
 */

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
			case '*' : Nucleotides[i] = '*';
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
