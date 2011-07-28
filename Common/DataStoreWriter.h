#ifndef _DATASTOREWRITER_H
#define _DATASTOREWRITER_H

#include <string>
#include "DataStore.h"

using namespace std;

class DataStoreWriter
{
public:
	DataStoreWriter();
	virtual ~DataStoreWriter();

public:
	bool Open(const string &fileName, const string &mode = "wb");
	bool Close();
	bool Write(const DataStore &store);

protected:
	FILE *out;
};
#endif
