#ifndef _TIMERS_H
#define _TIMERS_H

#include <map>
#include <sys/time.h>

using namespace std;

class Timers
{
public:
	Timers();
	virtual ~Timers();

public:
	int AddTimer();
	bool RemoveTimer(int id);
	timeval GetTimer(int id);
	double Elapsed(int id);

private:
	map<int, timeval> timers;
	int count;
};

#endif
