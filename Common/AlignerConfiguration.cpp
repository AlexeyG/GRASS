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
	SuffixArrayExactCommand = "bwa aln -n 0 -t %i -f %s %s %s >& /dev/null";
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
