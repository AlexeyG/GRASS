#ifndef _LOCATION_H
#define _LOCATION_H
#include <string>

using namespace std;

class Location
{
public:
	Location(int id, int coord);

public:
	string ToString() const;

public:
	int Id;
	int Coord;
};
#endif
