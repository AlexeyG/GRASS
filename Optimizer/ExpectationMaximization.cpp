#include "ExpectationMaximization.h"

ExpectationMaximization::ExpectationMaximization()
{
	status = Clean;
}

ExpectationMaximization::~ExpectationMaximization()
{
}

bool ExpectationMaximization::Formulate(const DataStore &store, const vector<double> &coord)
{
	//
}

bool ExpectationMaximization::Formulate(const DataStore &store)
{
	//
}

bool ExpectationMaximization::Solve()
{
	//
}

SolverStatus ExpectationMaximization::GetStatus() const
{
	return status;
}

double ExpectationMaximization::GetObjective() const
{
}

double ExpectationMaximization::GetDistanceSlack(int i) const
{
}

double ExpectationMaximization::GetOrderSlack(int i) const
{
}

int ExpectationMaximization::GetSlackCount() const
{
}
