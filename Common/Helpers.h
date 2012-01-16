/*
 * Common : a collection of classes (re)used throughout the scaffolder implementation.
 * Copyright (C) 2011  Alexey Gritsenko
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 * 
 * 
 * Email: a.gritsenko@tudelft.nl
 * Mail: Delft University of Technology
 *       Faculty of Electrical Engineering, Mathematics, and Computer Science
 *       Department of Mediamatics
 *       P.O. Box 5031
 *       2600 GA, Delft, The Netherlands
 */


#ifndef _HELPERS_H
#define _HELPERS_H

#include "DataStore.h"
#include "Timers.h"
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

namespace Helpers
{
	const double Eps = 1e-5;
	const double Inf = 1e10;
	const double Pi = atan2(0.0, -1.0);
	const string TempFilePrefix("tmp.");
	extern Timers ElapsedTimers;

	double RandomUniform();
	double RandomNormal(double mean = 0.0, double std = 1.0);
	int RandomNormal(int mean = 0, int std = 1);
	bool IsNumber(const char *str);
	string ItoStr(int a);
	int ParseInt(const char *str, bool &success);
	string TempFile(string path = "");
	string RandomString(int len, const char *alpha = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	bool FileExists(const string &fileName);
	bool RemoveFile(const string &fileName);
	bool Execute(const string &cmd);

	void PrintDataStore(const DataStore &store);
	template<class T>
	void Swap(T &a, T &b)
	{
		T t(a);
		a = b;
		b = t;
	}
        
        string NextEntry(string &str);
	template <class T>
	T GetArgument(const string &str)
        {
            stringstream ss(str);
            T res;
            ss >> res;
            return res;
        }
}
#endif
