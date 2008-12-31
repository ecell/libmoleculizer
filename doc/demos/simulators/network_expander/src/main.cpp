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
processCommandLineArgs( int argc, char* argv[], std::string& theFileName, int& number);

void
displayHelpAndExitProgram();

bool getUninitializedSpecies( const mzr::moleculizer& moleculizerRef, std::string& speciesName);

void printAllSpeciesStreams(mzr::moleculizer& theMolzer);

void printStreamByName( mzr::moleculizer& refMolzer, const std::string& streamName);
void printStreamByTag( mzr::moleculizer& refMolzer, const std::string& streamName);

void printAllSpeciesByName(mzr::moleculizer& theMolzer);
void printAllSpeciesByTag(mzr::moleculizer& theMolzer);

int main(int argc, char* argv[])
{

  std::string fileName;
  int number = -1;

  processCommandLineArgs(argc, argv, fileName, number);

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
          std::cout << "Iteration " << iterNdx << std::endl;

          std::string name;
          if (getUninitializedSpecies( theMoleculizer, name))
          {
              theMoleculizer.incrementNetworkBySpeciesName( name );
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

  std::cout << "################################################" << '\n';

  printAllSpeciesStreams(theMoleculizer);

  std::cout << "################################################" << '\n';

  printAllSpeciesByName(theMoleculizer);

  std::cout << "################################################" << '\n';

  printAllSpeciesByTag(theMoleculizer);

  std::cout << "################################################" << '\n';
  
  return 0;
  
}

void printAllSpeciesStreams(mzr::moleculizer& refMolzer)
{
    std::vector<std::string> theStreams;
    refMolzer.getSpeciesStreams( theStreams );
    std::cout << " There are " << theStreams.size() << " streams." << std::endl;

    BOOST_FOREACH( const std::string& name, theStreams)
    {
        printStreamByName( refMolzer, name);
        std::cout << "################################################" << '\n';
        printStreamByTag( refMolzer, name);
    }
}

void printStreamByName( mzr::moleculizer& refMolzer, const std::string& streamName)
{
    std::cout << "'" << streamName << "'\n[\n" ;

    std::vector<const mzr::mzrSpecies*> theSpecies;

    refMolzer.getSpeciesInSpeciesStream( streamName, theSpecies);

    BOOST_FOREACH( const mzr::mzrSpecies* ptrSpecies, theSpecies)
    {
        std::cout << streamName << "@@" << ptrSpecies->getName() << '\n';
    }

    std::cout << ']' << std::endl;
}

void printStreamByTag( mzr::moleculizer& refMolzer, const std::string& streamName)
{
    std::cout << "'" << streamName << "'\n[\n" ;

    std::vector<const mzr::mzrSpecies*> theSpecies;

    refMolzer.getSpeciesInSpeciesStream( streamName, theSpecies);

    BOOST_FOREACH( const mzr::mzrSpecies* ptrSpecies, theSpecies)
    {
        std::cout << streamName << "@@" << ptrSpecies->getTag() << '\n';
    }

    std::cout << ']' << std::endl;
}

bool getUninitializedSpecies( const mzr::moleculizer& moleculizerRef, std::string& speciesName)
{
    BOOST_FOREACH( const mzr::moleculizer::SpeciesCatalog::value_type& thePair, moleculizerRef.getSpeciesCatalog() )
    {
        if(!thePair.second->hasNotified())
        {
            speciesName = *thePair.first;
            return true;
        }
    }

    return false;
}


void processCommandLineArgs( int argc, char* argv[], std::string& mzrFile, int& number)
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
    }

    if ( !file )
    {
      std::cerr << "Error, a file must be specified with an -f <FILE> parameter." << std::endl;
      exit( 1 );
    }

}

void printAllSpeciesByName(mzr::moleculizer& theMolzer)
{
    BOOST_FOREACH( const mzr::moleculizer::SpeciesCatalog::value_type& vt, theMolzer.theSpeciesListCatalog)
        {
            std::cout << "ALL@@" << *vt.first << std::endl;
        }
}

void printAllSpeciesByTag(mzr::moleculizer& theMolzer)
{
  BOOST_FOREACH( const mzr::moleculizer::SpeciesCatalog::value_type& vt, theMolzer.theSpeciesListCatalog)
        {
            std::cout << "ALL@@" << vt.second->getTag() << std::endl;
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
