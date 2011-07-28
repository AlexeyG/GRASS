#ifndef _DATASTOREREADER_H
#define _DATASTOREREADER_H

#include <string>
#include <iostream>
#include "DataStore.h"

using namespace std;

class DataStoreReader
{
public:
	DataStoreReader() {};
	virtual ~DataStoreReader();

public:
	bool Open(const string &fileName);
	bool Close();
	bool Read(DataStore &store);

private:
	bool readHeader(int &nContigs, int &nGroups, int &nLinks);
	bool readContigs(int nContigs, DataStore &store);
	bool readGroups(int nGroups, DataStore &store);
	bool readLinks(int nLinks, DataStore &store);
	bool readContig(Contig &contig, int &id);
	bool readGroup(LinkGroup &group, int &id);
	bool readLink(int &groupID, ContigLink &link);

private:
	static string nextEntry(string &str);
	template <class T>
	static T getArgument(const string &str);

protected:
	fstream in;
};
#endif
