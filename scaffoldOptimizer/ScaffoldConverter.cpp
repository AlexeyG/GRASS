#include "ScaffoldConverter.h"
#include <string>
#include "Helpers.h"
#include "MinMax.h"

using namespace std;

FastASequence ScaffoldConverter::ToFasta(const DataStore &store, const Scaffold &scaffold)
{
	string name, sequence;
	int count = scaffold.ContigCount();
	int end = 0;
	for (int i = 0; i < count; i++)
	{
		ScaffoldContig contig = scaffold[i];
		FastASequence contigSeq = store[contig.Id].Sequence;
		int contigLen = contigSeq.Nucleotides.length();
		if (!contig.T)
		{
			cout << i << " Forward: " << contig.X << " - " << end << " = " << contig.X - end << endl;
			string spacer(max(contig.X - end, 0), 'N');
			sequence = sequence + spacer + contigSeq.Nucleotides;
			name = (!name.empty() ? name + "|+" : "+") + Helpers::ItoStr(contig.Id);
			end += spacer.length() + contigLen;
		}
		else
		{
			cout << i << " Reverse: " << contig.X << " - " << contigLen << " + 1 -" << end << " = " << contig.X - contigLen + 1 - end << endl;
			contigSeq.ReverseCompelement();
			string spacer(max(contig.X - contigLen + 1 - end, 0), 'N');
			sequence = sequence + spacer + contigSeq.Nucleotides;
			name = (!name.empty() ? name + "|-" : "-") + Helpers::ItoStr(contig.Id);
			end += spacer.length() + contigLen;
		}
	}

	return FastASequence(sequence, name);
}

vector<FastASequence> ScaffoldConverter::ToFasta(const DataStore &store, const vector<Scaffold> &scaffold)
{
	int count = scaffold.size();
	vector<FastASequence> seq(count);
	for (int i = 0; i < count; i++)
		seq[i] = ToFasta(store, scaffold[i]);

	return seq;
}
