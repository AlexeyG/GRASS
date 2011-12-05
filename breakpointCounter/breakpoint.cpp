/* 
 * File:   breakpoint.cpp
 * Author: alexeyg
 *
 * Created on December 5, 2011, 7:42 PM
 */

#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>
#include <boost/shared_array.hpp>
#include "Configuration.h"
#include "Sequence.h"
#include "Reader.h"
#include "Aligner.h"

using namespace std;

typedef vector<FastASequence> Sequences;

Configuration config;
auto_ptr<Sequences> scaffolds;
auto_ptr<Sequences> references;


bool readContigs(const string &fileName, Sequences &contigs)
{
    FastAReader reader;
    bool result = reader.Open(fileName) && reader.Read(contigs) > 0;
    reader.Close();
    return result;
}

bool alignScaffolds(const string &referenceFileName, const string &scaffoldsFilename)
{
    MummerAligner aligner(referenceFileName, scaffoldsFilename, config.MummerConfig);
    if (aligner.Align())
        cout << "Success: " << aligner.OutputFileName << endl;
    else
        cout << "Failure: " << aligner.OutputFileName << endl;
    return false;
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(NULL));
    scaffolds = auto_ptr<Sequences>(new Sequences());
    references = auto_ptr<Sequences>(new Sequences());
    if (config.ProcessCommandLine(argc, argv))
    {
        if (!readContigs(config.ScaffoldFileName, *scaffolds))
        {
            cerr << "[-] Unable to read scaffolds (" << config.ScaffoldFileName << ")." << endl;
            return -2;
        }
        cerr << "[+] Read scaffolds (" << config.ScaffoldFileName << ")." << endl;
        if (!readContigs(config.ReferenceFileName, *references))
        {
            cerr << "[-] Unable to read references (" << config.ReferenceFileName << ")." << endl;
            return -3;
        }
        cerr << "[+] Read references (" << config.ReferenceFileName << ")." << endl;
        if (!alignScaffolds(config.ReferenceFileName, config.ScaffoldFileName))
            return -4;
        cerr << "[+] Aligned scaffolds to reference (" << config.ScaffoldFileName << " -> " << config.ReferenceFileName << ")." << endl;
        
        /*if (!calculateDepth(*coverage, *contigs, *depth))
        {
            cerr << "[-] Unable to calculate coverage depth." << endl;
            return -4;
        }*/
        //
        return 0;
    }
    cerr << config.LastError;
    return -1;
}
