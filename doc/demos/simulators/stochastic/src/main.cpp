/////////////////////////////////////////////////////////////////////////////
// libMoleculizer - a species and reaction generator for reaction networks
// Copyright (C) 2008  Molecular Sciences Institute
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// This file was authored by Nathan Addy <addy@molsci.org> 
// Contact information:
//   Nathan Addy, Scientific Programmer     Email: addy@molsci.org
//   The Molecular Sciences Institute
//   2168 Shattuck Ave.
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include "utl/arg.hh"
#include "demostochasticsimulator.hpp"

using namespace std;

void 
processCommandLineArgs(int argc, char* argv[], std::string& theFileName, std::string& modelFile);
void displayHelpAndExitProgram();

int main(int argc, char* argv[])
{
    std::string modelfile;
    std::string rulesfile;
    processCommandLineArgs(argc, argv, rulesfile, modelfile);

    SimpleStochasticSimulator theSim( rulesfile, modelfile);

    std::cout << "\n\n#######################################################\nBeginning simulation\n\n";

    for(unsigned int ii = 0; ii != 500; ++ii)
    {
        theSim.singleStep();
    }

    std::cout << "\n\n#######################################################\nRun Completed.\n\nFinal State:\n";
    theSim.printState();
    return 0;

}

void processCommandLineArgs(int argc, char* argv[], std::string& rulesFile, std::string& modelFile)
{

    bool rulesFound(false);
    bool modelFound(false);
    // Skip the command name 
    argc--;
    argv++;

    if (argc == 0)
    {
        displayHelpAndExitProgram();
    }

    while(0 < argc)
    {
        std::string arg(*argv);
        argv++;
        argc--;

        if (arg == "--help")
        {
            displayHelpAndExitProgram();
        }
        
        if(arg == "--rules")
        {
            rulesFile = utl::mustGetArg(argc, argv);
            rulesFound = true;
        }

        if(arg == "--model")
        {
            modelFile = utl::mustGetArg(argc, argv);
            modelFound = true;
        }
    }

    if (!modelFound || !rulesFound)
    {
        std::cerr << "Error, both a --rules and a --model parameter must be provided." << endl;
        exit(1);
    }
    
}

void displayHelpAndExitProgram()
{

    cout << "Usage: stoch_demo --model modelfile --rules rulesfile" << endl;

    cout << "This is a demonstration program that demonstrates how libmoleculizer can be used as " << endl;
    cout << "a component in a Gillespie-like simulator." << endl;

    cout << "This program should have come with associated documentation.  Please read it for more details." << endl;
    cout << "\tNathan Addy <addy@molsci.org>\n\tSeptember 18, 2008." << endl;
    
    exit(0);
}
