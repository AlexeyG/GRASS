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
