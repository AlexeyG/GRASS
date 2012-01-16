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
