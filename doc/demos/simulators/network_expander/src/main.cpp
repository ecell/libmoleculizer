//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Nathan Addy, Scientific Programmer, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//


#include <string>
#include <iostream>
#include "utl/arg.hh"
#include "mzr/moleculizer.hh"
using namespace std;

void
processCommandLineArgs( int argc, char* argv[], std::string& theFileName);

void
displayHelpAndExitProgram();


int main(int argc, char* argv[])
{

  std::string fileName;

  processCommandLineArgs(argc, argv, fileName);

  mzr::moleculizer theMoleculizer;
  theMoleculizer.attachFileName( fileName );

  theMoleculizer.generateCompleteNetwork();

  std::cout << "There are " << theMoleculizer.getTotalNumberReactions() 
	    << " reactions and " << theMoleculizer.getTotalNumberSpecies() << " species in the nework." << std::endl;
  
  return 0;
  
}


void processCommandLineArgs( int argc, char* argv[], std::string& mzrFile)
{

    bool file( false );

    // Skip the command name
    argc--;
    argv++;

    if ( argc == 0 )
    {
        displayHelpAndExitProgram();
    }

    while ( 0 < argc )
    {
        std::string arg( *argv );
        argv++;
        argc--;

        if ( arg == "--help" )
        {
            displayHelpAndExitProgram();
        }

        if ( arg == "-f" || arg == "--file")
        {
	  mzrFile = utl::mustGetArg( argc, argv );
	  file = true;
        }
    }

    if ( !file )
    {
      std::cerr << "Error, a file must be specified with an -f <FILE> parameter." << std::endl;
      exit( 1 );
    }

}

void displayHelpAndExitProgram()
{

    cout << "Usage: network_expander -f <FILE>" << endl;

    cout << "This is a demonstration program that demonstrates how libmoleculizer can be used as " << endl;
    cout << "a component for expanding whole reaction networks and displaying some information." << endl;

    cout << "on them.";

    cout << "Libmoleculizer should have come with associated documentation.  Please read it for more details." << endl;
    cout << "\tNathan Addy <addy@molsci.org>\n\tSeptember 18, 2008." << endl;

    exit( 0 );
}
