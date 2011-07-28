#include "DataStoreWriter.h"
#include <cstdlib>
#include <cstdio>

DataStoreWriter::DataStoreWriter()
{
	out = NULL;
}

DataStoreWriter::~DataStoreWriter()
{
	if (out != NULL)
	{
		fclose(out);
		out = NULL;
	}
}

bool DataStoreWriter::Open(const string &fileName, const string &mode)
{
	if (out != NULL)
		return false;
	out = fopen(fileName.c_str(), mode.c_str());
	return out != NULL;
}

bool DataStoreWriter::Close()
{
	fclose(out);
	out = NULL;
	return true;
}

bool DataStoreWriter::Write(const DataStore &store)
{
	if (out == NULL)
		return false;
	int nContigs = store.ContigCount;
	int nGroups = store.GroupCount;
	int nLink = store.LinkCount;
	fprintf(out, "%i\t%i\t%i\n", nContigs, nGroups, nLink);
	for (int i = 0; i < nContigs; i++)
	{
		fprintf(out, "%i\t%s\n", store[i].GetID(), store[i].Sequence.Comment.c_str());
		fprintf(out, "%s\n", store[i].Sequence.Nucleotides.c_str());
	}
	for (int i = 0; i < nGroups; i++)
	{
		const LinkGroup &group = store.GetGroup(i);
		fprintf(out, "%i\t%s\t%s\n", group.GetID(), group.Name.c_str(), group.Description.c_str());
	}
	for (DataStore::LinkMap::const_iterator it = store.Begin(); it != store.End(); it++)
		fprintf(out, "%i\t%i\t%i\t%i\t%i\t%lf\t%lf\t%i\t%lf\t%s\n", it->second.GetGroupID(), it->first.first, it->first.second, (it->second.EqualOrientation ? 1 : 0), (it->second.ForwardOrder ? 1 : 0), it->second.Mean, it->second.Std, (it->second.Ambiguous ? 1 : 0), it->second.Weight, it->second.Comment.c_str());
	return true;
}
