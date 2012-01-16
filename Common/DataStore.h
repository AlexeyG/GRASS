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


#ifndef _DATASTORE_H
#define _DATASTORE_H

#include "Sequence.h"
#include <vector>
#include <string>
#include <map>
#include <set>

using namespace std;

class Contig
{
public:
	Contig(const FastASequence &seq = FastASequence()) : Sequence(seq), id(0) {};

public:
	int GetID() const;

public:
	FastASequence Sequence;

private:
	int id;

	friend class DataStore;
};

class ContigLink
{
public:
	ContigLink(int first = -1, int second = -1, double mean = 0, double std = 0, bool equalOrientation = false, bool forwardOrder = false, double weight = 0, const string &comment = string());

public:
	bool operator< (const ContigLink &other);
	int GetGroupID() const;

public:
	int First, Second;
	double Mean, Std;
	bool EqualOrientation;
	bool ForwardOrder;
	double Weight;
	bool Ambiguous;
	string Comment;

private:
	int groupId;

	friend class DataStore;
};

class LinkGroup
{
public:
	LinkGroup(const string &name, const string &description = (char *)"") : Name(name), Description(description), id(0) {};

public:
	int GetID() const;

public:
	string Name;
	string Description;

private:
	int id;

	friend class DataStore;
};

class DataStore
{
public:
	DataStore() : ContigCount(0), LinkCount(0), GroupCount(0) {};
	typedef multimap<pair<int,int>,ContigLink> LinkMap;
	typedef pair<LinkMap::const_iterator, LinkMap::const_iterator> LinkRange;

public:
	const Contig &operator[] (int i) const;
	const LinkRange operator() (int i, int j) const;
	LinkMap::const_iterator Begin() const;
	LinkMap::const_iterator End() const;
	LinkMap::iterator Begin();
	LinkMap::iterator End();
	const LinkGroup &GetGroup(int id) const;
        vector<FastASequence> GetContigs() const;
	int AddContig(const Contig &contig);
	int AddGroup(const LinkGroup &group);
	LinkMap::const_iterator AddLink(int groupId, const ContigLink &link);
	bool ReadContigs(const string &fileName);
	void Sort();
	void Bundle(bool sortLinks, bool perGroup, bool joinAmbiguous, double distance = 3);
	void Extract(const vector<int> &what, DataStore &store, vector<int> &transBack);
	int RemoveAmbiguous();
	int Erode(double weight);
        int IsolateContigs(const vector<int> &ids);

private:
	void bundleLinks(vector<ContigLink> &l, bool perGroup, bool joinAmbiguous, double distance);
	void performBundle(vector<ContigLink> &l, double distance);
	static bool linkComparer(const ContigLink &a, const ContigLink &b);
	static bool linkComparerAmbiguous(const ContigLink &a, const ContigLink &b);
	static bool linkComparerGroup(const ContigLink &a, const ContigLink &b);
	static bool linkComparerAmbiguousGroup(const ContigLink &a, const ContigLink &b);
	static bool selectGroup(vector<ContigLink> &l, bool perGroup, bool joinAmbiguous, int &s, vector<ContigLink> &selection);

public:
	int ContigCount;
	int LinkCount;
	int GroupCount;

private:
	vector<Contig> contigs;
	vector<LinkGroup> groups;
	LinkMap links;
};
#endif
