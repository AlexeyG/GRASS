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
	ContigLink(int first = -1, int second = -1, double mean = 0, double std = 0, bool equalOrientation = false, bool forwardOrder = false, double weight = 0);

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
	const LinkGroup &GetGroup(int id) const;
	int AddContig(const Contig &contig);
	int AddGroup(const LinkGroup &group);
	LinkMap::const_iterator AddLink(int groupId, const ContigLink &link);
	bool ReadContigs(const string &fileName);
	void Sort();
	void Bundle(bool sortLinks, bool perGroup, bool joinAmbiguous, double distance = 3);
	void Extract(const vector<int> &what, DataStore &store);

private:
	void bundleLinks(vector<ContigLink> &l, bool perGroup, bool joinAmbiguous, double distance);
	void performBundle(vector<ContigLink> &l, double distance);
	static bool linkComparer(const ContigLink &a, const ContigLink &b);
	static bool linkComparerGroup(const ContigLink &a, const ContigLink &b);
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
