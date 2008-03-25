/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008 Nathan Addy
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
/////////////////////////////////////////////////////////////////////////////

#include "moleculizer.hh"
#include "mzrHelperFunctions.hh"
#include "debug.hh" 

#include "utl/utlHelper.hh"
#include "utl/string.hh"
#include <algorithm>
#include <fstream>
#include <pause.hpp>

namespace mzr
{
  void 
  moleculizer::RunInteractiveDebugMode()
  {
    mzr::InteractiveModeManager<moleculizer> theInteractiveMode(this);

    theInteractiveMode.addFunction("showNumberSpecies", &moleculizer::DEBUG_showTotalNumberSpecies);
    theInteractiveMode.addFunction("showNumberReactions", &moleculizer::DEBUG_showTotalNumberReactions);
    theInteractiveMode.addFunction("showNumberDeltaSpecies", &moleculizer::DEBUG_showNumberDeltaSpecies);
    theInteractiveMode.addFunction("showNumberDeltaReactions", &moleculizer::DEBUG_showNumberDeltaReactions);
    theInteractiveMode.addFunction("showAllSpecies", &moleculizer::DEBUG_showAllSpecies);
    theInteractiveMode.addFunction("showAllReactions", &moleculizer::DEBUG_showAllReactions);
    theInteractiveMode.addFunction("showNewlyCreatedSpecies", &moleculizer::DEBUG_showDeltaSpecies);
    theInteractiveMode.addFunction("showNewlyCreatedReactions", &moleculizer::DEBUG_showDeltaReactions);
    theInteractiveMode.addFunction("showSingleLiveSpecies", &moleculizer::DEBUG_showRandomLiveSpecies);
    theInteractiveMode.addFunction("showNumberLiveSpecies", &moleculizer::DEBUG_showNumberLiveSpecies);
    theInteractiveMode.addFunction("incrementSpecies", &moleculizer::DEBUG_incrementSpecies);
    theInteractiveMode.addFunction("incrementRandomSpecies", &moleculizer::DEBUG_incrementRandomSpecies);
    theInteractiveMode.addFunction("calculate reaction between species", &moleculizer::DEBUG_showReactionsBetweenSpecies);

    theInteractiveMode.runInteractiveMode();
  }

  void moleculizer::RunProfileMode(unsigned int Num_Iterations, bool verbose)
  {
    for( unsigned int i = 0;
         i != Num_Iterations;
         ++i)
      {
        if (i % 10 == 0 )
          {
            cout << i / static_cast<float>(Num_Iterations)  * 100.0f << "% done..." << endl;
          }

        resetCurrentState();

        try
          {
            DEBUG_incrementNRandomLiveSpecies(1);
            DEBUG_showTotalNumberSpecies();
            DEBUG_showTotalNumberReactions();
          }
        catch(utl::xcpt x)
          {
            x.warn();
            cerr << "Unexpectedly could not find anything...Quitting" << endl;
            i = Num_Iterations;
            break;
          }
      }

    if (verbose )
      {
        DEBUG_showAllSpecies();
        DEBUG_showAllReactions();
      }
               
  }

  void
  moleculizer::DEBUG_showTotalNumberSpecies() 
  {
    cout << "There are " << getTotalNumberSpecies() << " species in the list." << endl;
  }

  void
  moleculizer::DEBUG_showTotalNumberReactions() 
  {
    cout << "There are " << getTotalNumberReactions() << " reactions in the list." << endl;
  }

  void 
  moleculizer::DEBUG_showNumberDeltaSpecies() 
  {
    cout << "There are " << getNumberDeltaSpecies() << " delta species." << endl;
  }

  void 
  moleculizer::DEBUG_showNumberDeltaReactions() 
  {
    cout << " There are " << getNumberDeltaReactions() << " delta reactions." << endl;
  }
    
  void
  moleculizer::DEBUG_showAllSpecies() 
  {
    unsigned int index = 0;

    for(SpeciesCatalogCIter iter = theSpeciesListCatalog.begin();
        iter != theSpeciesListCatalog.end();
        ++iter, ++index)
      {
        cout << *(iter->first) << endl;
      }
  }

  void 
  moleculizer::DEBUG_showAllReactions() 
  {
    std::for_each( theCompleteReactionList.begin(),
                   theCompleteReactionList.end(),
                   aux::printPtrWithName());
        
  }

  void
  moleculizer::DEBUG_showDeltaSpecies() 
  {
    std::for_each(theDeltaSpeciesList.begin(),
                  theDeltaSpeciesList.end(),
                  aux::printPtrWithName());
  }

  void moleculizer::DEBUG_showDeltaReactions() 
  {
    std::for_each(theDeltaReactionList.begin(),
                  theDeltaReactionList.end(),
                  aux::printPtrWithIndexedName());
  }

  void moleculizer::DEBUG_incrementSpecies()
  {
    resetCurrentState();

    std::string nameToIncrement;
    std::cin >> nameToIncrement;
        
    try
    {
        cout << "Incrementing '" << nameToIncrement << "'..." << endl;
        incrementNetworkBySpeciesName( nameToIncrement );
        cout << "Generation Results:" << endl;
        DEBUG_showNumberDeltaSpecies();
        DEBUG_showNumberDeltaReactions();
      }
    catch( utl::xcpt x)
      {
        cout << x.getMessage() << endl;
        cout << "Continuing...." << endl;
      }
  }

    
  std::string
  moleculizer::DEBUG_getRandomLiveSpecies() const 
  {

      SpeciesListCIter speciesListIter = theDeltaSpeciesList.begin();

      while(speciesListIter != theDeltaSpeciesList.end()) 
      {
          if (!(*speciesListIter)->hasNotified())
          {
              return (*speciesListIter)->getName();
          }

          ++speciesListIter;
      }

      SpeciesCatalogCIter speciesCatalogIter  = theSpeciesListCatalog.begin();
      while( speciesCatalogIter != theSpeciesListCatalog.end())
      {
          if ( !speciesCatalogIter->second->hasNotified() )
          {
              return *(speciesCatalogIter->first);
          }
          speciesCatalogIter++;
      }
             
      throw utl::xcpt("Error.  Reaction network has no \"live\" species.");
  }
  
  
  void
  moleculizer::DEBUG_incrementRandomSpecies()
  {
        std::string randomSpeciesName = DEBUG_getRandomLiveSpecies();
        cout << "Expanding '" << randomSpeciesName << "'..." << endl;
        incrementNetworkBySpeciesName( randomSpeciesName );
   
        cout << "Following expansion:" << endl;
        DEBUG_showNumberDeltaSpecies();
        DEBUG_showNumberDeltaReactions();
        DEBUG_showTotalNumberSpecies();
        DEBUG_showTotalNumberReactions();
  }
  
  void 
  moleculizer::DEBUG_incrementNRandomLiveSpecies(unsigned int numIters)
  {
    for (unsigned int index = 0;
         index != numIters;
         ++index)
      {
            DEBUG_incrementRandomSpecies();
      }
  }

  void moleculizer::DEBUG_showNumberLiveSpecies()
  {
    unsigned int number = 0;
    for( SpeciesCatalogCIter iter = theSpeciesListCatalog.begin();
         iter != theSpeciesListCatalog.end();
         ++iter)
      {
        if (!iter->second->hasNotified()) ++number;
      }
        
    cout << "There are " << number <<  " live species." << endl;
        
        
  }

  void moleculizer::DEBUG_showRandomLiveSpecies()
  {
    cout << DEBUG_getRandomLiveSpecies() << endl;
  }


  void moleculizer::DEBUG_showReactionsBetweenSpecies() 
  {

    std::string speciesOne, speciesTwo;
    cout << "-- Please enter two species names --\nSpecies 1:\t";
    cin >> speciesOne;
    cout << "\nSpecies 2:\t";
    cin >> speciesTwo;
    cout << endl;


    try
      {
        mzr::mzrSpecies* ptrSpeciesOne = getSpeciesFromName( speciesOne );
        mzr::mzrSpecies* ptrSpeciesTwo = getSpeciesFromName( speciesTwo );
      }
     catch(mzr::illegalSpeciesNameXcpt xcpt)
       {
       }
    catch(...)
      {
      }

  }



}




