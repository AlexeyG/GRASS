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


#include "AlignerConfiguration.h"

// Construtor with default configuration parameter settings.
BWAConfiguration::BWAConfiguration()
{
	TmpPath = "/tmp";
	NumberOfThreads = 8;
	MaximumHits = 1000;
	ExactMatch = false;
	IndexCommand = "bwa index -a is -p %s %s >& /dev/null";
	SuffixArrayCommand = "bwa aln -t %i -f %s %s %s >& /dev/null";
	SuffixArrayExactCommand = "bwa aln -t %i -f %s %s %s >& /dev/null";
	AlignSingleEndCommand = "bwa samse -f %s -n %i %s %s %s >& /dev/null";
}

// Construtor with default configuration parameter settings.
NovoAlignConfiguration::NovoAlignConfiguration()
{
	TmpPath = "/tmp";
	IndexCommand = "novoindex -m %s %s >& /dev/null";
	AlignSingleEndCommand = "(novoalign -d %s -f %s -o SAM -r All > %s) >& /dev/null";
}

// Construtor with default configuration parameter settings.
SAMToolsConfiguration::SAMToolsConfiguration()
{
	TmpPath = "/tmp";
	ConvertCommand = "samtools view -b -S -o %s %s >& /dev/null";
}

// Construction with default configuration parameter settings
MummerConfiguration::MummerConfiguration()
{
    TmpPath = "/tmp";
    NucmerCommand = "nucmer -p %s %s %s >& /dev/null";
    DeltaFilterCommand = "(delta-filter -q %s > %s) >& /dev/null";
    ShowCoordsCommand = "(show-coords -q -T -d -H %s > %s) >& /dev/null";
}

// Construction with default configuration parameter settings
MummerTilerConfiguration::MummerTilerConfiguration()
{
    TmpPath = "/tmp";
    NucmerCommand = "nucmer -p %s %s %s >& /dev/null";
    ShowTilingCommand = "(show-tiling %s > %s) >& /dev/null";
}
