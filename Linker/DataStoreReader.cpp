#include "DataStoreReader.h"
#include <sstream>

DataStoreReader::~DataStoreReader()
{
	if (in.is_open())
		in.close();
}

bool DataStoreReader::Open(const string &fileName)
{
	if (in.is_open())
		return false;
	in.open(fileName, ios::in);
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
	string contigsStr = nextEntry(line);
	string groupsStr = nextEntry(line);
	string linksStr = nextEntry(line);
	if (contigsStr.length() == 0 || groupsStr.length() == 0 || linksStr.length() == 0)
		return false;
	nContigs = getArgument<int>(contigsStr);
	nGroups = getArgument<int>(groupsStr);
	nLinks = getArgument<int>(linksStr);
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
	int groupIdD;
	for (int i = 0; i < nGroups; i++)
	{
		if (in.eof() || !readGroup(group, groupIdD))
			return false;
		if (store.AddGroup(group) != groupIdD)
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
}

bool DataStoreReader::readContig(Contig &contig, int &id)
{
	string line, seq;
	getline(in, line);
	if (in.eof()) return false;
	getline(in, seq);
	string idStr = nextEntry(line);
	string commentStr = nextEntry(line);
	seq = nextEntry(seq);
	if (idStr.length() == 0 || commentStr.length() == 0 || seq.length() == 0)
		return false;
	id = getArgument<int>(idStr);
	contig = Contig(FastASequence(seq, commentStr));
	return true;
}

bool DataStoreReader::readGroup(LinkGroup &group, int &id)
{
	string line;
	getline(in, line);
	string idStr = nextEntry(line);
	string nameStr = nextEntry(line);
	string descriptionStr = nextEntry(line);
	if (idStr.length() == 0 || nameStr.length() == 0)
		return false;
	id = getArgument<int>(idStr);
	group = LinkGroup(nameStr, descriptionStr);
	return true;
}

bool DataStoreReader::readLink(int &groupID, ContigLink &link)
{
	string line;
	getline(in, line);
	string groupIdStr = nextEntry(line);
	string firstStr = nextEntry(line);
	string secondStr = nextEntry(line);
	string orientationStr = nextEntry(line);
	string orderStr = nextEntry(line);
	string meanStr = nextEntry(line);
	string stdStr = nextEntry(line);
	string ambiguousStr = nextEntry(line);
	string weightStr = nextEntry(line);
	if (firstStr.length() == 0 || secondStr.length() == 0 || orientationStr.length() == 0 || orderStr.length() == 0 || meanStr.length() == 0 || stdStr.length() == 0 || ambiguousStr.length() == 0 || weightStr.length() == 0)
		return false;
	groupID = getArgument<int>(groupIdStr);
	int first = getArgument<int>(firstStr);
	int second = getArgument<int>(secondStr);
	bool equalOrientation = getArgument<int>(orientationStr) == 1;
	bool forwardOrder = getArgument<int>(orderStr) == 1;
	double mean = getArgument<double>(meanStr);
	double std = getArgument<double>(stdStr);
	bool ambiguous = getArgument<int>(ambiguousStr) == 1;
	double weight = getArgument<double>(weightStr);
	if (groupID < 0 || first < 0 || second < 0 || weight <= 0)
		return false;
	link = ContigLink(first, second, mean, std, equalOrientation, forwardOrder, weight);
	link.Ambiguous = ambiguous;
	return true;
}

string DataStoreReader::nextEntry(string &str)
{
	int len = str.length();
	if (len == 0)
		return string();
	int pos = 0;
	while (pos < len && str[pos] != '\t' && str[pos] != '\n')
		pos++;
	string result = str.substr(0, pos);
	str = str.substr(pos + 1, len - pos - 1);
	return result;
}

template <class T>
T getArgument(const string &str)
{
	stringstream ss(str);
	T res;
	ss >> res;
	return res;
}
