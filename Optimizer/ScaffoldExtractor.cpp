#include "ScaffoldExtractor.h"
#include "DPGraph.h"
#include "Helpers.h"

ScaffoldContig::ScaffoldContig(int id, bool t, double x)
	: Id(id), T(t), X(x)
{
}

bool ScaffoldContig::operator< (const ScaffoldContig &b) const
{
	return X < b.X;
}

Scaffold::Scaffold()
{
}

void Scaffold::AddContig(int id, bool t, double x)
{
	AddContig(ScaffoldContig(id, t, x));
}

void Scaffold::AddContig(const ScaffoldContig &contig)
{
	contigs.push_back(contig);
}

void Scaffold::Sort()
{
	sort(contigs.begin(), contigs.end());
}

const ScaffoldContig &Scaffold::operator[] (int i) const
{
	return contigs[i];
}

int Scaffold::ContigCount() const
{
	return contigs.size();
}

vector<Scaffold> ScaffoldExtractor::Extract(const DataStore &store, bool single)
{
	vector<Scaffold> ans;
	if (single)
		extractSingleScaffold(store, ans);
	else
		extractOrientedScaffold(store, ans);
	return ans;
}

vector<Scaffold> ScaffoldExtractor::Extract(const IterativeSolver &solver)
{
	vector<Scaffold> ans;
	int m = solver.GetSlackCount();
	vector< vector<int> > components;
	vector<bool> slacks(m, true);
	for (int i = 0; i < m; i++)
		if (solver.GetDistanceSlack(i) >= FixedMIQPSolver::DesiredDistanceSlackMax && solver.GetOrderSlack(i) >= FixedMIQPSolver::DesiredOrderSlackMax)
			slacks[i] = false;
	DPGraph graph(solver.GetStore(), solver.T, slacks);
	graph.FindConnectedComponents(components);
	for (int i = 0; i < (int)components.size(); i++)
	{
		int size = components[i].size();
		Scaffold scaffold;
		for (int j = 0; j < size; j++)
			scaffold.AddContig(components[i][j], solver.T[components[i][j]], solver.X[components[i][j]]);
		scaffold.Sort();
		ans.push_back(scaffold);
	}
	return ans;
}

vector<Scaffold> ScaffoldExtractor::Extract(const DataStore &store, const GASolver &solver)
{
	vector<Scaffold> ans;
	vector< vector<int> > components;
	DPGraph graph(store, solver.T);
	graph.FindConnectedComponents(components);
	for (int i = 0; i < (int)components.size(); i++)
	{
		int size = components[i].size();
		Scaffold scaffold;
		for (int j = 0; j < size; j++)
			scaffold.AddContig(components[i][j], solver.T[components[i][j]], solver.X[components[i][j]]);
		scaffold.Sort();
		ans.push_back(scaffold);
	}
	return ans;
}

vector<Scaffold> ScaffoldExtractor::Extract(const DataStore &store, const FixedMIQPSolver &solver)
{
	vector<Scaffold> ans;
	int m = solver.GetSlackCount();
	vector< vector<int> > components;
	vector<bool> slacks(m, true);
	for (int i = 0; i < m; i++)
		if (solver.GetDistanceSlack(i) >= FixedMIQPSolver::DesiredDistanceSlackMax && solver.GetOrderSlack(i) >= FixedMIQPSolver::DesiredOrderSlackMax)
			slacks[i] = false;
	DPGraph graph(store, solver.T, slacks);
	graph.FindConnectedComponents(components);
	for (int i = 0; i < (int)components.size(); i++)
	{
		int size = components[i].size();
		Scaffold scaffold;
		for (int j = 0; j < size; j++)
			scaffold.AddContig(components[i][j], solver.T[components[i][j]], solver.X[components[i][j]]);
		scaffold.Sort();
		ans.push_back(scaffold);
	}
	return ans;
}

void ScaffoldExtractor::extractSingleScaffold(const DataStore &store, vector<Scaffold> &ans)
{
	vector< vector<int> > components;
	DPGraph graph(store);
	graph.FindConnectedComponents(components);
	for (int i = 0; i < (int)components.size(); i++)
	{
		int size = components[i].size();
		Scaffold scaffold;
		for (int j = 0; j < size; j++)
		{
			bool isForward;
			double position;
			getOrientation(store[components[i][j]], isForward, position);
			scaffold.AddContig(components[i][j], isForward, position);
		}
		scaffold.Sort();
		ans.push_back(scaffold);
	}
}

void ScaffoldExtractor::extractOrientedScaffold(const DataStore &store, vector<Scaffold> &ans)
{
	int n = store.ContigCount;
	vector<bool> t(n, false);
	vector<double> x(n, 0);
	for (int i = 0; i < n; i++)
	{
		bool val;
		getOrientation(store[i], val, x[i]);
		t[i] = val;
	}
	vector< vector<int> > components;
	DPGraph graph(store, t);
	graph.FindConnectedComponents(components);
	for (int i = 0; i < (int)components.size(); i++)
	{
		int size = components[i].size();
		Scaffold scaffold;
		for (int j = 0; j < size; j++)
			scaffold.AddContig(components[i][j], t[components[i][j]], x[components[i][j]]);
		scaffold.Sort();
		ans.push_back(scaffold);
	}
}

bool ScaffoldExtractor::getOrientation(const Contig &contig, bool &orientation, double &position)
{
	string name = contig.Sequence.Name();
	if (name[0] == '+')
		orientation = false;
	else if (name[0] == '-')
		orientation = true;
	else
		return false;
	position = 0;
	for (int i = 1; i < (int)name.length(); i++)
		if (name[i] == '|')
			break;
		else
			position = position * 10 + name[i] - '0';
	return true;
}
