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

#include "XATag.h"
#include "Helpers.h"

XATag::XATag(const BamAlignment &alg)
{
	RefID = alg.RefID;
	Position = alg.Position;
	IsReverseStrand = alg.IsReverseStrand();
}

XATag::XATag(const string &str, const BamReader &reader)
{
        if (str.empty())
                throw new exception();
        int n = str.length();
        int i = 0, s = 0;
        while (i < n - 1 && str[i] != ',') i++;
        string name = str.substr(0, i);
        s = ++i; while (i < n - 1 && str[i] != ',') i++;
        string position = str.substr(s, i - s);
        if (position.empty())
                throw new exception();
        if (position[0] != '+' && position[0] != '-')
                throw new exception();
        RefID = reader.GetReferenceID(name);
		bool tmp;
        Position = Helpers::ParseInt(position.c_str() + 1, tmp);
        IsReverseStrand = position[0] == '-';
}
