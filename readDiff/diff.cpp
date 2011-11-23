#include "Configuration.h"
#include "Reader.h"
#include "Writer.h"
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cstdlib>

using namespace std;

Configuration config;
vector<FastQSequence> A, B;

bool isNameLess(const FastQSequence &a, const FastQSequence &b)
{
	return strcmp(a.Comment.c_str(), b.Comment.c_str()) < 0;
}

bool readSet(const string &fileName, vector<FastQSequence> &reads)
{
	FastQReader reader;
	bool result = reader.Open(fileName) && reader.Read(reads) > 0;
	reader.Close();
	return result;
}

bool calculateDifference(const vector<FastQSequence> &A, const vector<FastQSequence> &B, const string &fileName)
{
	FastQWriter writer;
	bool result = writer.Open(fileName);
	int i = 0, j = 0;
	int aSize = A.size(), bSize = B.size();
	while (result && (i < aSize || j < bSize))
	{
		if (j >= bSize)
			result = result && writer.Write(A[i++]);
		else if (i >= aSize)
			j++;
		else
		{
			int cmp = strcmp(A[i].Comment.c_str(), B[j].Comment.c_str());
			if (cmp == 0)
				i++, j++;
			else if (cmp < 0)
				result = result && writer.Write(A[i++]);
			else
				j++;
		}
	}
	writer.Close();
	return result;
}

int main(int argc, char *argv[])
{
	if (config.ProcessCommandLine(argc, argv))
	{
		if (!readSet(config.AFileName, A))
		{
			cerr << "[-] Unable to read first input file: " << config.AFileName << endl;
			return -2;
		}
		cerr << "[+] Successfully read first input file." << endl;
		if (!readSet(config.BFileName, B))
		{
			cerr << "[-] Unable to read first input file: " << config.BFileName << endl;
			return -3;
		}
		sort(A.begin(), A.end(), isNameLess);
		sort(B.begin(), B.end(), isNameLess);
		cerr << "[+] Successfully read second input file." << endl;
		cerr << "[i] Calculating difference..." << endl;
		if (!calculateDifference(A, B, config.CFileName))
		{
			cerr << "[-] Unable to calculate the difference." << endl;
			return -4;
		}
		cerr << "[+] Successfully wrote difference to file: " << config.CFileName << endl;
		return 0;
	}
	cerr << config.LastError;
	return -1;
}
