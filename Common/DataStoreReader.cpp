#include "DataStoreReader.h"
#include "Helpers.h"
#include <sstream>

#include <iostream>

using namespace std;

DataStoreReader::~DataStoreReader()
{
	if (in.is_open())
		in.close();
}

bool DataStoreReader::Open(const string &fileName)
{
	if (in.is_open())
		return false;
	in.open(fileName.c_str(), ios::in);
	return in.is_open();
}

bool DataStoreReader::Close()
{
	in.close();
	return true;
}

bool DataStoreReader::Read(DataStore &store)
{
	int nContigs, nGroups, nLinks;
	if (!in.is_open())
		return false;
	if (!readHeader(nContigs, nGroups, nLinks))
		return false;
	if (!readContigs(nContigs, store))
		return false;
	if (!readGroups(nGroups, store))
		return false;
	if (!readLinks(nLinks, store))
		return false;
	return true;
}

bool DataStoreReader::readHeader(int &nContigs, int &nGroups, int &nLinks)
{
	string line;
	getline(in, line);
	string contigsStr = Helpers::NextEntry(line);
	string groupsStr = Helpers::NextEntry(line);
	string linksStr = Helpers::NextEntry(line);
	if (contigsStr.length() == 0 || groupsStr.length() == 0 || linksStr.length() == 0)
		return false;
	nContigs = Helpers::GetArgument<int>(contigsStr);
	nGroups = Helpers::GetArgument<int>(groupsStr);
	nLinks = Helpers::GetArgument<int>(linksStr);
	if (nContigs <= 0 || nGroups <= 0 || nLinks < 0)
		return false;
	return true;
}

bool DataStoreReader::readContigs(int nContigs, DataStore &store)
{
	Contig contig;
	int contigID;
	for (int i = 0; i < nContigs; i++)
	{
		if (in.eof() || !readContig(contig, contigID))
			return false;
		if (store.AddContig(contig) != contigID)
			return false;
	}
	return true;
}

bool DataStoreReader::readGroups(int nGroups, DataStore &store)
{
	LinkGroup group("");
	int groupId;
	for (int i = 0; i < nGroups; i++)
	{
		if (in.eof() || !readGroup(group, groupId))
			return false;
		if (store.AddGroup(group) != groupId)
			return false;
	}
	return true;
}

bool DataStoreReader::readLinks(int nLinks, DataStore &store)
{
	int groupID;
	ContigLink link;
	for (int i = 0; i < nLinks; i++)
	{
		if (in.eof() || !readLink(groupID, link))
			return false;
		store.AddLink(groupID, link);
	}
	return true;
}

bool DataStoreReader::readContig(Contig &contig, int &id)
{
	string line, seq;
	getline(in, line);
	if (in.eof()) return false;
	getline(in, seq);
	string idStr = Helpers::NextEntry(line);
	string commentStr = Helpers::NextEntry(line);
	seq = Helpers::NextEntry(seq);
	if (idStr.length() == 0 || commentStr.length() == 0 || seq.length() == 0)
		return false;
	id = Helpers::GetArgument<int>(idStr);
	contig = Contig(FastASequence(seq, commentStr));
	return true;
}

bool DataStoreReader::readGroup(LinkGroup &group, int &id)
{
	string line;
	getline(in, line);
	string idStr = Helpers::NextEntry(line);
	string nameStr = Helpers::NextEntry(line);
	string descriptionStr = Helpers::NextEntry(line);
	if (idStr.length() == 0 || nameStr.length() == 0)
		return false;
	id = Helpers::GetArgument<int>(idStr);
	group = LinkGroup(nameStr, descriptionStr);
	return true;
}

bool DataStoreReader::readLink(int &groupID, ContigLink &link)
{
	string line;
	getline(in, line);
	string groupIdStr = Helpers::NextEntry(line);
	string firstStr = Helpers::NextEntry(line);
	string secondStr = Helpers::NextEntry(line);
	string orientationStr = Helpers::NextEntry(line);
	string orderStr = Helpers::NextEntry(line);
	string meanStr = Helpers::NextEntry(line);
	string stdStr = Helpers::NextEntry(line);
	string ambiguousStr = Helpers::NextEntry(line);
	string weightStr = Helpers::NextEntry(line);
	string commentStr = Helpers::NextEntry(line);
	if (firstStr.length() == 0 || secondStr.length() == 0 || orientationStr.length() == 0 || orderStr.length() == 0 || meanStr.length() == 0 || stdStr.length() == 0 || ambiguousStr.length() == 0 || weightStr.length() == 0)
		return false;
	groupID = Helpers::GetArgument<int>(groupIdStr);
	int first = Helpers::GetArgument<int>(firstStr);
	int second = Helpers::GetArgument<int>(secondStr);
	bool equalOrientation = Helpers::GetArgument<int>(orientationStr) == 1;
	bool forwardOrder = Helpers::GetArgument<int>(orderStr) == 1;
	double mean = Helpers::GetArgument<double>(meanStr);
	double std = Helpers::GetArgument<double>(stdStr);
	bool ambiguous = Helpers::GetArgument<int>(ambiguousStr) == 1;
	double weight = Helpers::GetArgument<double>(weightStr);
	if (groupID < 0 || first < 0 || second < 0 || weight <= 0)
		return false;
	link = ContigLink(first, second, mean, std, equalOrientation, forwardOrder, weight, commentStr);
	link.Ambiguous = ambiguous;
	return true;
}
