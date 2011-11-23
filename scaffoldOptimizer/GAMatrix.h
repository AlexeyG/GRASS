#ifndef _GAMATRIX_H
#define _GAMATRIX_H
#include <vector>

using namespace std;

class GAMatrix
{
public:
	typedef vector< vector<double> > Matrix;

public:
	GAMatrix(int n = 0);

public:
	void CalculatePositions();

public:
	vector<double> &operator[] (int i);
	const vector<double> &operator[] (int i) const;

private:
	int size;
	Matrix matrix;

public:
	vector< vector<int> > Pos;
};
#endif
