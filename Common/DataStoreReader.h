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


#ifndef _DATASTOREREADER_H
#define _DATASTOREREADER_H

#include <string>
#include <fstream>
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

protected:
	fstream in;
};
#endif
