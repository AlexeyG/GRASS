#include "Helpers.h"
#include "Globals.h"
#include <string>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include<iostream>

using namespace std;

double Helpers::RandomUniform()
{
	double a = rand();
	return a / RAND_MAX;
}

double Helpers::RandomNormal(double mean, double variance)
{
	double u = RandomUniform();
	double v = RandomUniform();
	return mean + variance * sqrt(-2 * log(u)) * sin(2 * atan2(0.0, -1.0) * v);
}

int Helpers::RandomNormal(int mean, int variance)
{
	return (int)Helpers::RandomNormal((double)mean, (double)variance);
}

bool Helpers::IsNumber(char *str)
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

int Helpers::ParseInt(char *str, bool &success)
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
