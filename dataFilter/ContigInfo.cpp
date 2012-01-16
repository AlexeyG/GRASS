/*
 * dataFilter : a support tool used in debuging to filter paired reads aligning
 * to gaps between contigs.
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
