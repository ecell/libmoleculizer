//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2009 The Molecular Sciences Institute.
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
processCommandLineArgs( int argc, char* argv[], std::string& theFileName, int& number, int& printOutput, std::string& theOutputFileName);

void
displayHelpAndExitProgram();

bool getUninitializedSpecies( const mzr::moleculizer& moleculizerRef, std::string& speciesName);

void printAllSpeciesStreams(mzr::moleculizer& theMolzer);

void printStreamByName( mzr::moleculizer& refMolzer, const std::string& streamName);
void printStreamByTag( mzr::moleculizer& refMolzer, const std::string& streamName);

void printAllSpeciesByName(mzr::moleculizer& theMolzer);
void printAllSpeciesByID(mzr::moleculizer& theMolzer);

int main(int argc, char* argv[])
{

  std::string fileName;
  std::string outputFileName("");
  int number = -1;
  int printOutput = 1;

  processCommandLineArgs(argc, argv, fileName, number, printOutput, outputFileName);

  mzr::moleculizer theMoleculizer;
  theMoleculizer.attachFileName( fileName );

  if (number < 0 )
  {
      theMoleculizer.generateCompleteNetwork();
  }
  else
  {
      for ( int iterNdx = 0; iterNdx != number; ++iterNdx)
      {
          std::cout << "Iteration " << iterNdx << "/" << number << std::endl;

          std::string name;
          if (getUninitializedSpecies( theMoleculizer, name))
          {
              theMoleculizer.incrementNetworkBySpeciesTag( name );
          }
          else
          {
              std::cout << "All species incremented." << std::endl;
              std::cout << "Breaking on iteration " << iterNdx << std::endl;
              break;
          }
      }
  }

  std::cout << "################################################" << '\n';
  std::cout << "After " << number << " iterations," << '\n';
  std::cout << "There are " 
            << theMoleculizer.getTotalNumberReactions() << " reactions and " 
            << theMoleculizer.getTotalNumberSpecies() << " species in " 
            << theMoleculizer.getNumberOfPlexFamilies() << " families in the nework." << std::endl;

  if (printOutput)
  {

      std::cout << "################################################" << '\n';

      printAllSpeciesStreams(theMoleculizer);

      std::cout << "################################################" << '\n';

      printAllSpeciesByName(theMoleculizer);

      std::cout << "################################################" << '\n';

      printAllSpeciesByID(theMoleculizer);

      std::cout << "################################################" << '\n';
  }

  if(outputFileName != "")
  {
      theMoleculizer.writeOutputFile(outputFileName, true);
  }

  return 0;
  
}

void printAllSpeciesStreams(mzr::moleculizer& refMolzer)
{
    std::vector<std::string> theStreams;
    refMolzer.getSpeciesStreams( theStreams );
    std::cout << " There are " << theStreams.size() << " streams." << std::endl;

    for( std::vector<std::string>::const_iterator strIter = theStreams.begin();
         strIter != theStreams.end();
         ++strIter)
    {
        printStreamByName( refMolzer, *strIter);
        std::cout << "################################################" << '\n';
        printStreamByTag( refMolzer, *strIter);
    }
}

void printStreamByName( mzr::moleculizer& refMolzer, const std::string& streamName)
{
    std::cout << "'" << streamName << "'\n[\n" ;

    std::vector<const mzr::mzrSpecies*> theSpecies;

    refMolzer.getSpeciesInSpeciesStream( streamName, theSpecies);

    for( std::vector<const mzr::mzrSpecies*>::const_iterator specIter = theSpecies.begin();
         specIter != theSpecies.end();
         ++specIter)
    {
        std::cout << streamName << "@@" << (*specIter)->getName() << '\n';
    }

    std::cout << ']' << std::endl;
}

void printStreamByTag( mzr::moleculizer& refMolzer, const std::string& streamName)
{
    std::cout << "'" << streamName << "'\n[\n" ;

    std::vector<const mzr::mzrSpecies*> theSpecies;

    refMolzer.getSpeciesInSpeciesStream( streamName, theSpecies);

    for( std::vector<const mzr::mzrSpecies*>::const_iterator specIter = theSpecies.begin();
         specIter != theSpecies.end();
         ++specIter)
    {
        std::cout << streamName << "@@" << (*specIter)->getTag() << '\n';
    }

    std::cout << ']' << std::endl;
}

bool getUninitializedSpecies( const mzr::moleculizer& moleculizerRef, std::string& speciesName)
{
    for( mzr::moleculizer::SpeciesCatalog::const_iterator specIter = moleculizerRef.getSpeciesCatalog().begin();
         specIter != moleculizerRef.getSpeciesCatalog().end();
         ++specIter)
    {
        if(!specIter->second->hasNotified())
        {
            speciesName = *(specIter->first);
            return true;
        }
    }

    return false;
}


void processCommandLineArgs( int argc, char* argv[], std::string& mzrFile, int& number, int& print, std::string& outputFileName)
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
        if( arg == "-n" )
        {
            std::string numAsString = utl::mustGetArg( argc, argv);
            number = utl::argMustBeNNInt( numAsString );
        }
        if( arg == "-v" )
        {
            print = 1;
        }
        if(arg == "-q")
        {
            print = 0;
        }
        if(arg == "-o")
        {
            outputFileName = utl::mustGetArg( argc, argv );
        }
    }

    if ( !file )
    {
      std::cerr << "Error, a file must be specified with an -f <FILE> parameter." << std::endl;
      exit( 1 );
    }

}

void printAllSpeciesByName(mzr::moleculizer& theMolzer)
{
    for( mzr::moleculizer::SpeciesCatalog::const_iterator specIter = theMolzer.theSpeciesListCatalog.begin();
         specIter != theMolzer.theSpeciesListCatalog.end();
         ++specIter)
    {
        std::cout << "ALL@@" << *specIter->first << std::endl;
    }
}

void printAllSpeciesByID(mzr::moleculizer& theMolzer)
{
    for( mzr::moleculizer::SpeciesCatalog::const_iterator specIter = theMolzer.theSpeciesListCatalog.begin();
         specIter != theMolzer.theSpeciesListCatalog.end();
         ++specIter)
    {
        std::cout << "ALL@@" << theMolzer.convertSpeciesTagToSpeciesID( *specIter->first ) << std::endl;
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
