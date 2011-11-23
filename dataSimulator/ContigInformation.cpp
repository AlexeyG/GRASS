#include "ContigInformation.h"
#include <string>
#include <sstream>

using namespace std;

string ContigInformation::FormatName() const
{
	stringstream str;
	str << (ReverseOrientation ? "-" : "+") << Position << "|";
	return str.str();
}
