#include "Location.h"
#include <sstream>

using namespace std;

Location::Location(int id, int coord)
	: Id(id), Coord(coord)
{
}

string Location::ToString() const
{
	stringstream sstr;
	sstr << Id << "." << Coord;
	return sstr.str();
}
