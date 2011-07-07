#ifndef _HELPERS_H
#define _HELPERS_H
namespace Helpers
{
	double RandomUniform();
	double RandomNormal(double mean = 0.0, double variance = 1.0);
	int RandomNormal(int mean = 0, int variance = 1);
	bool IsNumber(char *str);
	int ParseInt(char *str, bool &success);
}
#endif
