#include "SolverConfiguration.h"

SolverConfiguration::SolverConfiguration()
{
	UseOpportunisticSearch = true;
	UseObjectiveHeuristic = true;
	SuppressOutput = false;
	VerboseOutput = false;
	Threads = 0;
	TimeLimit = 0;
	LPThreads = 0;
	LPTimeLimit = 30;
	LPAttempts = 5;
	GARestarts = 0;
}
