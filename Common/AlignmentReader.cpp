#include "AlignmentReader.h"
#include <string>
#include <cstdlib>

using namespace std;

AlignmentReader::AlignmentReader()
{
	buffer.Name.clear();
}

bool AlignmentReader::Open(const string &fileName)
{
	return reader.Open(fileName);
}

bool AlignmentReader::Close()
{
	reader.Close();
	return true;
}

bool AlignmentReader::IsOpen() const
{
	return reader.IsOpen();
}

bool AlignmentReader::GetNextAlignmentGroup(vector<BamAlignment> &alg)
{
	alg.clear();
	if (!buffer.Name.empty())
	{
		alg.push_back(buffer);
		buffer.Name.clear();
	}
	bool read;
	while ((read = reader.GetNextAlignment(buffer)))
	{
		if (!alg.empty() && (alg.end() - 1)->Name != buffer.Name)
			break;
		alg.push_back(buffer);
	}
	if (!read) buffer.Name.clear();
	return !alg.empty();
}

bool AlignmentReader::GetNextAlignmentGroup(BamAlignment &alignment, vector<XATag> &tags)
{
	tags.clear();
	vector<BamAlignment> alg;
	if (!GetNextAlignmentGroup(alg))
		return false;
	alignment = alg[0];
	for (vector<BamAlignment>::iterator it = alg.begin(); it != alg.end(); it++)
	{
		if (!it->IsMapped())
			continue;
		tags.push_back(XATag(*it));
		string xaStr;
		it->GetTag("XA", xaStr);
		parseXATag(xaStr, tags);
	}
	return true;
}

void AlignmentReader::parseXATag(const string &str, vector<XATag> &tags)
{
	int offset = 0;
	int found = (int)string::npos;
	while ((found = str.find(';', offset)) != (int)string::npos)
	{
		tags.push_back(parseXAEntry(str.substr(offset, found - offset)));
		offset = found + 1;
	}
}

XATag AlignmentReader::parseXAEntry(const string &str)
{
	int n = str.length();
	int i = 0;
	while (i < n - 1 && str[i] != ',') i++;
	string name = str.substr(0, i);
	int s = ++i;
	while (i < n - 1 && str[i] != ',') i++;
	string position = str.substr(s, i - s);
	return XATag(atoi(position.c_str() + 1), reader.GetReferenceID(name), position[0] == '-');
}
