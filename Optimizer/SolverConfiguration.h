#ifndef _SOLVERCONFIGURATION_H
#define _SOLVERCONFIGURATION_H
struct SolverConfiguration
{
public:
	SolverConfiguration();
public:
	bool UseOpportunisticSearch;
	bool UseObjectiveHeuristic;
	bool SuppressOutput;
	int VerboseOutput;
	int TimeLimit;
	int Threads;
	int LPThreads;
	int LPTimeLimit;
	int LPAttempts;
	int GATimeLimit;
	int GARestarts;
};
#endif
