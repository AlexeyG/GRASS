#include "Helpers.h"
#include "Globals.h"
#include <string>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

#include<iostream>

using namespace std;

double Helpers::RandomUniform()
{
	double a = rand();
	return a / RAND_MAX;
}

double Helpers::RandomNormal(double mean, double std)
{
	double u = RandomUniform();
	double v = RandomUniform();
	return mean + std * sqrt(-2 * log(u)) * sin(2 * Helpers::Pi * v);
}

int Helpers::RandomNormal(int mean, int std)
{
	return (int)Helpers::RandomNormal((double)mean, (double)std);
}

bool Helpers::IsNumber(const char *str)
{
	int len = strlen(str);
	if (len == 0)
		return false;
	
	int i = 0;
	if (str[0] == '-')
		i++;

	if (len == i)
		return false;

	for (; i < len; i++)
		if (!isdigit(str[i]))
			return false;
	
	return true;
}

string Helpers::ItoStr(int a)
{
	char buf[MaxLine];
	sprintf(buf, "%i", a);
	return string(buf);
}

int Helpers::ParseInt(const char *str, bool &success)
{
	int len = strlen(str);
	success = Helpers::IsNumber(str);
	if (!success)
		return 0;

	int i = 0;
	int sgn = +1;
	if (str[0] == '-')
	{
		sgn = -1;
		i++;
	}
	if (str[0] == '+')
		i++;
	int num = 0;
	while (i < len && str[i] == '0')
		i++;
	while (i < len)
	{
		num = num * 10 + str[i] - '0';
		i++;
	}
	return sgn * num;
}

string Helpers::RandomString(int len, const char *alpha)
{
	int n = strlen(alpha);
	string res;
	res.resize(len);
	while (len > 0)
	{
		len--;
		res[len] = alpha[rand() % n];
	}
	return res;
}

string Helpers::TempFile(string path)
{
	string fileName;
	if (path.length() > 0 && *path.end() != '/')
		path += "/";
	do
	{
		fileName = path + Helpers::TempFilePrefix + Helpers::RandomString(6);
	} while (Helpers::FileExists(fileName));
	return fileName;
}

bool Helpers::FileExists(const string &fileName)
{
	struct stat s;
	//cout << "Checking " << fileName << " : " << stat(fileName.c_str(), &s) << endl;
	return stat(fileName.c_str(), &s) == 0;
}

bool Helpers::RemoveFile(const string &fileName)
{
	return remove(fileName.c_str()) == 0;
}

bool Helpers::Execute(const string &cmd)
{
	return system(cmd.c_str()) == 0;
}

void Helpers::PrintDataStore(const DataStore &store)
{
	int n = store.ContigCount;
	int totalCnt = 0;
	double totalW = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double w = 0;
			int cnt = 0;
			const DataStore::LinkRange range = store(i, j);
			for (DataStore::LinkMap::const_iterator it = range.first; it != range.second; it++)
				w += it->second.Weight, cnt++;
			printf("%3i (%6.2lf)\t", cnt, w);
			totalCnt += cnt, totalW += w;
		}
		printf("\n");
	}
	printf("Sum: %i (%6.2lf)\n", totalCnt, totalW);
}

string Helpers::NextEntry(string &str)
{
	int len = str.length();
	if (len == 0)
		return string();
	int pos = 0;
	while (pos < len && str[pos] != '\t' && str[pos] != '\n')
		pos++;
	string result = str.substr(0, pos);
	if (len - pos - 1 <= 0 || pos + 1 >= len)
		str = string();
	else
		str = str.substr(pos + 1, len - pos - 1);
	return result;
}

Timers Helpers::ElapsedTimers;
