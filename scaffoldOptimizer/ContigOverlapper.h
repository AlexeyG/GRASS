/*
 * scaffoldOptimizer : solves the MIQP optimization and produces linear scaffold
 * sequences.
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


/* 
 * File:   ContigOverlapper.h
 * Author: alexeyg
 *
 * Created on December 11, 2011, 3:39 PM
 */

#ifndef _CONTIGOVERLAPPER_H
#define	_CONTIGOVERLAPPER_H

#include "OverlapperConfiguration.h"
#include <string>

using namespace std;

class ContigOverlapper
{
public:
    static int FindBestOverlap(const std::string &left, const std::string &right, int distance, const OverlapperConfiguration &config, double &score, std::string &consensus);
    static int GetOverlapSequences(int offset, const std::string &left, const std::string &right, std::string &leftSequence, std::string &rightSequence);
    
private:
    static int GetAlignmentScore(const std::string &left, const std::string &right);
    static std::string GetConsensus(const std::string &left, const std::string &right);
};

#endif	/* _CONTIGOVERLAPPER_H */
