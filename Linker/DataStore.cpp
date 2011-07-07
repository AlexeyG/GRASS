#include "DataStore.h"
#include "Reader.h"
#include "Helpers.h"
#include <string>
#include <utility>
#include <algorithm>
#include <set>

#include <iostream>

using namespace std;

int Contig::GetID() const
{
	return id;
}

ContigLink::ContigLink(int first, int second, double mean, double std, bool equalOrientation, bool forwardOrder, double weight)
	: First(first), Second(second), Mean(mean), Std(std), EqualOrientation(equalOrientation), ForwardOrder(forwardOrder), Weight(weight), Ambiguous(false), groupId(0)
{
}

bool ContigLink::operator< (const ContigLink &other)
{
	if (First < other.First)
		return true;
	if (First > other.First)
		return false;
	if (Second < other.Second)
		return true;
	if (Second > other.Second)
		return true;
	return false;
}

int ContigLink::GetGroupID() const
{
	return groupId;
}

int LinkGroup::GetID() const
{
	return id;
}

const Contig &DataStore::operator[] (int i) const
{
	return contigs[i];
}

const DataStore::LinkRange DataStore::operator() (int i, int j) const
{
	pair<int,int> query(i, j);
	return links.equal_range(query);
}

const LinkGroup &DataStore::GetGroup(int id) const
{
	return groups[id];
}

int DataStore::AddContig(const Contig &contig)
{
	contigs.push_back(contig);
	int id = contigs[ContigCount].id = ContigCount++;
	return id;
}

int DataStore::AddGroup(const LinkGroup &group)
{
	groups.push_back(group);
	return groups[GroupCount].id = GroupCount++;
}

DataStore::LinkMap::const_iterator DataStore::AddLink(int groupId, const ContigLink &link)
{
	if (groupId >= GroupCount)
		throw exception();
	LinkCount++;
	fprintf(stderr, "Linked: d(%i,%i) = %8.2f; orientation: %8s; order: %7s.\n", link.First, link.Second, link.Mean, (link.EqualOrientation ? "equal" : "opposite"), (link.ForwardOrder ? "forward" : "reverse"));
	pair<int,int> query(link.First, link.Second);
	return links.insert(pair<pair<int,int>,ContigLink>(query, link));
}

bool DataStore::ReadContigs(const string &fileName)
{
	FastAReader reader;
	FastASequence seq;
	if (!reader.Open(fileName))
		return false;
	
	int read = 0;
	while (reader.Read(seq))
	{
		Contig contig(seq);
		AddContig(contig);
		read++;
	}
	reader.Close();
	return read != 0;
}

void DataStore::Sort()
{
	vector<ContigLink> vec;
	for (LinkMap::iterator it = links.begin(); it != links.end(); it++)
	{
		if (it->first.first > it->first.second)
		{
			/*ContigLink old = it->second;
			bool order = (old.EqualOrientation ? !old.ForwardOrder : old.ForwardOrder);
			ContigLink link(old.Second, old.First, old.Mean, old.Std, old.EqualOrientation, order, old.Weight);
			link.Ambiguous = old.Ambiguous;*/
			swap(it->second.First, it->second.Second);
			if (!it->second.EqualOrientation)
				it->second.ForwardOrder = !it->second.ForwardOrder;
		}
		vec.push_back(it->second);
	}
	links.clear();
	for (vector<ContigLink>::iterator it = vec.begin(); it != vec.end(); it++)
	{
		pair<int, int> pos(it->First, it->Second);
		links.insert(pair<pair<int,int>,ContigLink>(pos, *it));
	}
}

void DataStore::Bundle(bool sortLinks, bool perGroup, bool joinAmbiguous, double distance)
{
	if (sortLinks)
		Sort();

	vector<ContigLink> vec;
	vector< vector<ContigLink> > space;
	for (LinkMap::iterator it = links.begin(), start = links.begin(); it != links.end(); it++)
	{
		vec.clear();
		while (it != links.end() && it->first == start->first)
		{
			vec.push_back(it->second);
			it++;
		}
		start = it--;
		space.push_back(vec);
	}
	links.clear();
	for (vector< vector<ContigLink> >::iterator it = space.begin(); it != space.end(); it++)
		bundleLinks(*it, perGroup, joinAmbiguous, distance);
}

void DataStore::bundleLinks(vector<ContigLink> &l, bool perGroup, bool joinAmbiguous, double distance)
{
	if (perGroup)
		sort(l.begin(), l.end(), linkComparerGroup);
	else
		sort(l.begin(), l.end(), linkComparer);
	int s = 0;
	vector<ContigLink> sel;
	while (selectGroup(l, perGroup, joinAmbiguous, s, sel))
		while (sel.size() > 0)
			performBundle(sel, distance);
}

void DataStore::performBundle(vector<ContigLink> &l, double distance)
{
	vector<ContigLink> r;
	set<int> groupIds;
	int n = l.size(), m = (n % 2 == 1 ? n / 2 : n / 2 - 1);
	if (n == 0) return;
	double p = 0, q = 0, w = 0;
	bool ambiguous = false;
	for (int i = 0; i < n; i++)
		if (abs(l[i].Mean - l[m].Mean) < distance * l[m].Std)
		{
			p += l[i].Mean / (l[i].Std * l[i].Std);
			q += 1 / (l[i].Std * l[i].Std);
			w += l[i].Weight;
			ambiguous = ambiguous || l[i].Ambiguous;
			groupIds.insert(l[i].GetGroupID());
		}
		else r.push_back(l[i]);
	ContigLink link(l[m].First, l[m].Second, p / q, 1 / sqrt(q), l[m].EqualOrientation, l[m].ForwardOrder, w);
	link.Ambiguous = ambiguous;
	int groupId;
	if (groupIds.size() > 1)
	{
		string name;
		string description = "Group for bundle-links joining several groups.";
		for (set<int>::iterator it = groupIds.begin(); it != groupIds.end(); it++)
			name = name + "+" + Helpers::ItoStr(*it);
		groupId = AddGroup(LinkGroup(name, description));
	}
	else groupId = *groupIds.begin();
	AddLink(groupId, link);
	l.clear();
	l.insert(l.begin(), r.begin(), r.end());
}

bool DataStore::linkComparer(const ContigLink &a, const ContigLink &b)
{
	if (!a.EqualOrientation && b.EqualOrientation)
		return true;
	if (a.EqualOrientation && !b.EqualOrientation)
		return false;
	if (!a.ForwardOrder && b.ForwardOrder)
		return true;
	if (a.ForwardOrder && !b.ForwardOrder)
		return false;
	if (!a.Ambiguous && b.Ambiguous)
		return true;
	if (a.Ambiguous && !b.Ambiguous)
		return false;
	return a.Mean < b.Mean;
}

bool DataStore::linkComparerGroup(const ContigLink &a, const ContigLink &b)
{
	if (a.GetGroupID() < b.GetGroupID())
		return true;
	if (a.GetGroupID() > b.GetGroupID())
		return false;
	return linkComparer(a, b);
}

bool DataStore::selectGroup(vector<ContigLink> &l, bool perGroup, bool joinAmbiguous, int &s, vector<ContigLink> &selection)
{
	int i = s;
	int n = l.size();
	if (i >= n)
		return false;
	s++;
	selection.clear();
	selection.push_back(l[i]);
	while (s < n && ((!perGroup || l[i].GetGroupID() == l[s].GetGroupID()) && l[i].EqualOrientation == l[s].EqualOrientation && l[i].ForwardOrder == l[s].ForwardOrder && (joinAmbiguous || l[i].Ambiguous == l[s].Ambiguous)))
	{
		selection.push_back(l[s]);
		s++;
	}
	return true;
}
