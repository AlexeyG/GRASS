#include "Timers.h"

Timers::Timers()
	: count(0)
{
}

Timers::~Timers()
{
}

int Timers::AddTimer()
{
	timeval t;
	gettimeofday(&t, NULL);
	timers[count++] = t;
	return count - 1;
}

bool Timers::RemoveTimer(int id)
{
	map<int, timeval>::iterator it = timers.find(id);
	if (it == timers.end())
		return false;
	timers.erase(it);
	return true;
}

timeval Timers::GetTimer(int id)
{
	return timers[id];
}

double Timers::Elapsed(int id)
{
	timeval u = timers[id], v;
	gettimeofday(&v, NULL);
    double elapsedTime = (v.tv_sec - u.tv_sec) * 1000.0;
    elapsedTime += (v.tv_usec - u.tv_usec) / 1000.0;

	return elapsedTime;
}
