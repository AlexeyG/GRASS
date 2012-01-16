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


#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <string>
#include "api/BamAlignment.h"

using namespace std;
using namespace BamTools;

class Sequence
{
public:
	Sequence() {};
	Sequence(const string &sequence) : Nucleotides(sequence) {};
	Sequence(const BamAlignment &alg);

public:
	virtual void ReverseCompelement();

public:
	string Nucleotides;

protected:
	void complement();
};

class FastASequence : public Sequence
{
public:
	FastASequence() : Sequence() {};
	FastASequence(const string &sequence, const string &comment) : Sequence(sequence), Comment(comment) { };
	FastASequence(const BamAlignment &alg);

public:
	string Name() const;

public:
	string Comment;
};

class FastQSequence : public FastASequence
{
	public:
	FastQSequence() : FastASequence() {};
	FastQSequence(const string &sequence, const string &comment, const string &quality) : FastASequence(sequence, comment) { Quality = quality; };
	FastQSequence(const BamAlignment &alg);

public:
	virtual void ReverseCompelement();
public:
	string Quality;
};

#endif
