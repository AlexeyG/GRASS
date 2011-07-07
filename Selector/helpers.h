#ifndef _HELPERS_H
#define _HELPERS_H

#include <cmath>
#include <string>

using namespace std;

namespace Helpers
{
	const double Eps = 1e-10;
	const double Pi = atan2(0, -1);
	const string TempFilePrefix("tmp.");

	double RandomUniform();
	double RandomNormal(double mean = 0.0, double variance = 1.0);
	int RandomNormal(int mean = 0, int variance = 1);
	bool IsNumber(char *str);
	int ParseInt(char *str, bool &success);
	string TempFile(string path = "");
	string RandomString(int len, char *alpha = (char *)"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	bool FileExists(const string &fileName);
}
#endif
