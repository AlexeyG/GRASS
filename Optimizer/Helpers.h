#ifndef _HELPERS_H
#define _HELPERS_H

#include "DataStore.h"
#include "Timers.h"
#include <cmath>
#include <string>

using namespace std;

namespace Helpers
{
	const double Eps = 1e-5;
	const double Inf = 1e10;
	const double Pi = atan2(0.0, -1.0);
	const string TempFilePrefix("tmp.");
	extern Timers ElapsedTimers;

	double RandomUniform();
	double RandomNormal(double mean = 0.0, double variance = 1.0);
	int RandomNormal(int mean = 0, int variance = 1);
	bool IsNumber(char *str);
	string ItoStr(int a);
	int ParseInt(char *str, bool &success);
	string TempFile(string path = "");
	string RandomString(int len, const char *alpha = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	bool FileExists(const string &fileName);
	bool RemoveFile(const string &fileName);
	bool Execute(const string &cmd);

	void PrintDataStore(const DataStore &store);
	template<class T>
	void Swap(T &a, T &b)
	{
		T t(a);
		a = b;
		b = t;
	}
}
#endif
