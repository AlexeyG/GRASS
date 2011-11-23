#ifndef _GASOLVER_H
#define _GASOLVER_H
#include "Solver.h"
#include "DataStore.h"
#include "GAIndividual.h"
#include "GAMatrix.h"
#include <vector>

class GAIndividual;
class GASolver : public Solver
{
public:
	GASolver();
	virtual ~GASolver();

public:
	virtual bool Formulate(const DataStore &store);
	bool Formulate(const DataStore &store, const vector<double> &distanceSlack, const vector<double> &orderSlack);
	virtual bool Solve();
	virtual SolverStatus GetStatus() const;
	virtual double GetObjective() const;

public:
	void AddIndividual(const vector<bool> &t);

private:
	bool shouldTerminate();
	int generatePopulation(int from = 0);
	int localSearch(int from = 0);
	int crossover();
	void select();
	int restart(int from = 0);
	void formulateMatrix(const DataStore &store, const vector<double> &distanceSlack, const vector<double> &orderSlack);
	void selectInitialSolution();
	void updateSolution(const GAIndividual &ind);
	double getTime(double &lastIteration) const;

public:
	int SelectionSize;
	int LocalSearchM;
	int RestartGenerations;
	double CrossoverRate;

protected:
	int timerId;
	int iteration;
	int restartCount;
	int lastSuccess;
	double bestObjective;
	int populationSize;
	vector<GAIndividual> population;
public: // remove me!
	GAMatrix matrix;

public:
	GAIndividual InnovativeCrossover(const GAIndividual &p1, const GAIndividual &p2);
	void RandomizedKopt(GAIndividual &ind);
	void Mutate(GAIndividual &ind);

	//bool checkObjective(GAIndividual &ind);
	//bool checkGains(GAIndividual &ind);
};

#endif
