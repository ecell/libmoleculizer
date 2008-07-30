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
#include "nmr/nmrUnit.hh"
#include "mzr/unitsMgr.hh"

#include "utl/debug.hh"
#include "utl/utlHelper.hh"
#include "utl/utility.hh"
#include <algorithm>
#include <fstream>

namespace mzr
{
  void 
  moleculizer::RunInteractiveDebugMode()
  {
    mzr::InteractiveModeManager<moleculizer> theInteractiveMode(this);

    theInteractiveMode.addFunction("showLiveSpecies", &moleculizer::DEBUG_showLiveSpecies);
    theInteractiveMode.addFunction("showDeadSpecies", &moleculizer::DEBUG_showDeadSpecies);
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
    theInteractiveMode.addFunction("Change naming strategy", &moleculizer::DEBUG_changeNamingStrategy);
    theInteractiveMode.addFunction("Get species from name", &moleculizer::DEBUG_getSpeciesFromName);
    theInteractiveMode.addFunction("Clear all species and reactions", &moleculizer::DEBUG_clearAll);
    theInteractiveMode.addFunction("Find an unary reaction", &moleculizer::DEBUG_showUnaryReaction);
    theInteractiveMode.addFunction("Do random particle collision", &moleculizer::DEBUG_doRandomParticleCollisionInterval);

    theInteractiveMode.runInteractiveMode();
  }

  void moleculizer::DEBUG_doRandomParticleCollisionInterval()
  {
      std::cout << "DEBUG_doRandomParticleCollisionInterval not yet implemented." << endl;
      return;

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
  moleculizer::DEBUG_getRandomDeadSpeciesName() const
  {
      std::vector<mzrSpecies*> deadSpecies;
      BOOST_FOREACH( SpeciesCatalog::value_type vt, theSpeciesListCatalog)
      {
          if (vt.second->hasNotified())
          {
              deadSpecies.push_back( vt.second );
          }
      }

      if (deadSpecies.size() == 0) throw utl::xcpt("Error.  Reaction network has no \"live\" species.");
      
      int random_index = rand() % deadSpecies.size();

      mzrSpecies* pLiveSpecies = deadSpecies[random_index];
      std::string name = pLiveSpecies->getName();

      return name;
  }

  std::string
  moleculizer::DEBUG_getRandomLiveSpeciesName() const 
  {

      // Highly inefficient - this could be rewritten in a much
      // much better way by caching the liveSpecies. Each iteration
      // would begin by clearing the deltaSpecies list.
      // Next a random element would be chosen, removed from 
      // the vector, and incremented.  Then, any deltaSpecies would
      // be copied to the end of the vector.  
      // This would be FAR better...
      std::vector<mzrSpecies*> liveSpecies;
      BOOST_FOREACH( SpeciesCatalog::value_type vt, theSpeciesListCatalog)
      {
          if (!vt.second->hasNotified())
          {
              liveSpecies.push_back( vt.second );
          }
      }

      if (liveSpecies.size() == 0) throw utl::xcpt("Error.  Reaction network has no \"live\" species.");
      
      int random_index = rand() % liveSpecies.size();
      mzrSpecies* pLiveSpecies = liveSpecies[random_index];
      std::string name = pLiveSpecies->getName();
      return name;
  }
  
  
  void
  moleculizer::DEBUG_incrementRandomSpecies()
  {
      // Clear the delta lists....
      resetCurrentState();

        std::string randomSpeciesName = DEBUG_getRandomLiveSpeciesName();
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
    cout << DEBUG_getRandomLiveSpeciesName() << endl;
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
           mzr::mzrSpecies* pSpeciesOne = this->getSpeciesWithName( speciesOne );
           mzr::mzrSpecies* pSpeciesTwo = this->getSpeciesWithName( speciesTwo );

           std::vector<mzrReaction*> rxnVector;
           bool reactionBetweenSpecies = this->findReactionWithSubstrates(pSpeciesOne, pSpeciesTwo, rxnVector);

         if(reactionBetweenSpecies)
           {
               std::cout << "Found " << rxnVector.size() << " reactions between these" << std::endl;

               BOOST_FOREACH(mzrReaction* ptr, rxnVector)
               {
                   cout << ptr->getName() << endl;
               }
           }
         else
           {
               std::cout << "No reaction between" << std::endl;
           }

       }
     catch(mzr::illegalSpeciesNameXcpt xcpt)
     {
         xcpt.warn();
         std::cerr << "Continuing." << std::endl;
     }
  }



    void moleculizer::DEBUG_showLiveSpecies()
    {

      BOOST_FOREACH( SpeciesCatalog::value_type vt, theSpeciesListCatalog)
      {
          if (!vt.second->hasNotified())
          {
              cout << *vt.first << endl;
          }
      }
    }

    void moleculizer::DEBUG_showDeadSpecies()
    {
      BOOST_FOREACH( SpeciesCatalog::value_type vt, theSpeciesListCatalog)
      {
          if (vt.second->hasNotified())
          {
              cout << *vt.first << endl;
          }
      }
    }
    

  void moleculizer::DEBUG_showUnaryReaction() 
  {

    std::string speciesOne, speciesTwo;
    cout << "-- Please enter a species name --\nSpecies 1:\t";
    cin >> speciesOne;
    cout << endl;

    int count = 0;

    if ( theSpeciesListCatalog.find( &speciesOne) == theSpeciesListCatalog.end())
    {
        pUserUnits->pNmrUnit->constructSpeciesFromName( speciesOne );
        if (theSpeciesListCatalog.find( &speciesOne) == theSpeciesListCatalog.end() )
        {
            cerr << "error" << endl;
            exit(1);
        }
    }
    
    incrementNetworkBySpeciesName( speciesOne );
        
    std::vector<mzrReaction*> aVector;

    mzrSpecies* ptrSpecies = theSpeciesListCatalog[&speciesOne];
    findReactionWithSubstrates(ptrSpecies , aVector);

    std::sort(aVector.begin(),
              aVector.end());

    std::vector<mzrReaction*>::iterator newend = std::unique(aVector.begin(),
                                                             aVector.end());

    std::vector<mzrReaction*> bVector(aVector.begin(),
                                      newend);

    cout << "Found:\t" << aVector.size() << endl;

    BOOST_FOREACH(mzrReaction* ptr, bVector)
    {
        cout << ptr->getName() << endl;
    }

  }

  void moleculizer::DEBUG_changeNamingStrategy()
  {
    unsigned int choice = 0;

    while( choice < 1 || choice > 3)
      {
        cout << "1:\tDefault (NOT implemented yet)" << endl;
        cout << "2:\tInformative (NOT implemented yet)" << endl;
        cout << "3:\tCompact" << endl;
        cout << "Pick a naming scheme: ";
        cin >> choice;
      }

    std::string newNameEncoder("");

    if( choice == 1)
      {
        newNameEncoder = "basic-name-assembler";
      }
    else if (choice == 2)
      {
        newNameEncoder = "detailed-name-assembler";
      }
    else if (choice == 3)
      {
        newNameEncoder = "mangled-name-assembler";
      }
    
    pUserUnits->pNmrUnit->setDefaultNameEncoder( newNameEncoder );

    SpeciesCatalog newCatalog;
    
    for(SpeciesCatalog::const_iterator specIter = theSpeciesListCatalog.begin();
        specIter != theSpeciesListCatalog.end();
        ++specIter)
      {

        
        mzrSpecies* ptrMzrSpecies = specIter->second;
        delete specIter->first;

        std::string* newName = new std::string( ptrMzrSpecies->getName() );
        newCatalog.insert( std::make_pair(newName, ptrMzrSpecies));
      }
    
    theSpeciesListCatalog.clear();

    theSpeciesListCatalog.insert( newCatalog.begin(), newCatalog.end() );




  }

void 
moleculizer::DEBUG_doMultipleRandomParticleCollisions(unsigned int numCollisions)
{
      
}


void
moleculizer::DEBUG_getSpeciesFromName()
{
    cout << "Enter species name:\t" ;
    std::string name;
    cin >> name;

    mzr::mzrSpecies* pSpecies;
    try
      {
        pSpecies = getSpeciesWithName( name );
      }
    catch(IllegalNameXcpt xcpt)
      {
        xcpt.warn();
        return;
      }

    if(pSpecies)
      {

          std::string outputName = pSpecies->getName();
    
          cout << "You got back species with name:\t'" << pSpecies->getName() << "'"<< endl;

          if (name == outputName) 
          {
              std::cout << "they are the same." << std::endl;
          }
          else
          {
              std::cout << "They are NOT the same" << std::endl;
          }
      }
    else
      {
	std::cout << "That was an ILLEGAL name." << std::endl;
      }
}


    void
    moleculizer::DEBUG_isPresent() const
    {
        mzr::mzrSpecies* ptrMzr = (mzr::mzrSpecies*) 0xfffffff;
        BOOST_FOREACH( SpeciesCatalog::value_type type, theSpeciesListCatalog)
        {
            if (type.second == ptrMzr)
            {
                std::cout << "Present" << endl;
                return;
            }
        }
        
        std::cout << "Not present." << endl;
    }

    void 
    moleculizer::DEBUG_clearAll()
    {
        // First, we can delete all the species that have been created 
    
    }



}




