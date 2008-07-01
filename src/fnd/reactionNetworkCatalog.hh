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


#ifndef RXNNETWORKCATALOG_HH
#define RXNNETWORKCATALOG_HH


#include "fndExceptions.hh"
#include "rnHelperClasses.hh"
#include "utl/autoCatalog.hh"
#include "utl/macros.hh"
#include "fnd/basicReaction.hh"
#include "plex/mzrPlexFamily.hh"
#include "plex/mzrPlexSpecies.hh"
#include "utl/debug.hh"

#include <boost/foreach.hpp>
#include <vector>
#include <list>
#include <stack>
#include <sstream>
#include <iostream>
#include <string>

namespace fnd
{

  // I obviously need to add some kind of public API function for expanding by a speciesID.

  // Both speciesType and reactionType must have a public function
  // void expandReactionNetwork(). 
  // That is,  derived from fnd::reactionNetworkComponent.
  
  template <typename speciesT,
	    typename reactionT>
  class ReactionNetworkDescription
  {
    
  public:
    DECLARE_TYPE( speciesT, SpeciesType);
    DECLARE_TYPE( reactionT, ReactionType);

    typedef std::map< const std::string*, SpeciesTypePtr, aux::compareByPtrValue<std::string> > SpeciesCatalog;
    typedef typename SpeciesCatalog::iterator SpeciesCatalogIter;
    typedef typename SpeciesCatalog::const_iterator SpeciesCatalogCIter;

    typedef std::list<SpeciesTypePtr> SpeciesList;
    typedef std::list<ReactionTypePtr> ReactionList;

    typedef std::multimap<SpeciesTypePtr, ReactionTypePtr> ParticipatingSpeciesRxnMap;

    typedef std::list< std::pair<std::string*, SpeciesTypePtr> > SpeciesListCatalog;

    typedef typename SpeciesListCatalog::iterator SpeciesListCatalogIter;
    typedef typename SpeciesListCatalog::const_iterator SpeciesListCatalogCIter;
        
    typedef typename SpeciesList::iterator SpeciesListIter;
    typedef typename SpeciesList::const_iterator SpeciesListCIter;
    
    typedef typename ReactionList::iterator ReactionListIter;
    typedef typename ReactionList::const_iterator ReactionListCIter;

  public:

    SpeciesTypePtr
    findSpecies(const std::string& name) throw(fnd::NoSuchSpeciesXcpt)
    {
      SpeciesCatalogIter theIter = theSpeciesListCatalog.find(&name);
      if (theIter != theSpeciesListCatalog.end() )
	{
	  return theIter->second;
	}
      else
	{
	  throw fnd::NoSuchSpeciesXcpt(name);
	}
            
    }

    SpeciesTypeCptr
    findSpecies(const std::string& name) const throw(fnd::NoSuchSpeciesXcpt)
    {
      SpeciesCatalogCIter theIter = theSpeciesListCatalog.find(&name);
      
      if (theIter != theSpeciesListCatalog.end() )
	{
	  return theIter->second;
	}
      else
	{
	  throw fnd::NoSuchSpeciesXcpt(name);
	}
    }


    void
    findReactionWithSubstrates(const SpeciesTypePtr A, 
			       std::vector<ReactionTypePtr>& reactionVector)
    {

        reactionVector.clear();
        
        typename ParticipatingSpeciesRxnMap::const_iterator iter = singleSubstrateRxns.find(A);
        

        while( iter->first == A)

        {
            reactionVector.push_back(iter->second);
            ++iter;
        }

        return;
    }
      
      void
      findReactionWithSubstrates(const SpeciesTypePtr A, 
                                 const SpeciesTypePtr B,
                                 std::vector<ReactionTypePtr>& reactionVector) const
      {

          reactionVector.clear();
          const SpeciesTypePtr tmp;

          // This can be optimized by searching off the ptrSpececiesType which is a substrate
          // of fewer reactions.
          if (doubleSubstrateRxns.count(A) < doubleSubstrateRxns.count(B) )
          {
              for(typename ParticipatingSpeciesRxnMap::const_iterator iter = doubleSubstrateRxns.lower_bound(A);
                  iter != doubleSubstrateRxns.upper_bound(A);
                  ++iter)
              {
                  if ( iter->second->getReactants().find(B) != (iter->second->getReactants().end()) )
                  {
                      reactionVector.push_back( iter->second );
                  }
              }
                      
          }
          else
          {
              for(typename ParticipatingSpeciesRxnMap::const_iterator iter = doubleSubstrateRxns.lower_bound(B);
                  iter != doubleSubstrateRxns.upper_bound(B);
                  ++iter)
              {
                  if (iter->second->getReactants().find(A) != iter->second->getReactants().end())
                  {
                      reactionVector.push_back( iter->second);
                  }
              }
          } 
              

          return;

      }




    // The unary reaction case
    // This will ONLY return an unary reaction.  
    const ReactionTypePtr
    findReactionWithSubstrates(const SpeciesTypePtr A) const
    {

      int count = 0;
      std::string reactantName = A->getName();

      ReactionTypePtr thePtr = NULL;

      // TODO Finish writing this function.
      for(ReactionListCIter iter = theCompleteReactionList.begin();
	  iter != theCompleteReactionList.end();
	  ++iter)
	{
	  if ((*iter)->getReactants().size() != 1) continue;
	  else 
	    {
	      const SpeciesTypePtr B = ((*iter)->getReactants().begin()->first);
                        
	      if (reactantName == B->getName())
		{
		  count += 1;
		  thePtr = *iter;
		}
	    }
                

	}
      return thePtr;
    }

    const ReactionTypePtr
    findReactionWithSubstrates( SpeciesTypePtr A, SpeciesTypePtr B) const
    {
      for(ReactionListCIter iter = theCompleteReactionList.begin();
	  iter != theCompleteReactionList.end();
	  ++iter)
	{
	  ReactionTypePtr pRxn = *iter;

	  bool seenA = false;
	  bool seenB = false;
                      
	  for( typename std::map<SpeciesType*, int>::const_iterator jiter = pRxn->getReactants().begin();
	       jiter != pRxn->getReactants().end();
	       ++jiter)
	    {
	      if( jiter->first == A )
		{
		  seenA = true;
		}
	      else if (jiter->first == B)
		{
		  seenB = true;
		}
	    }

	  if ( seenA && seenB )
	    {
	      return pRxn;
	    }

	}

      // Scanned through everything and found nothing.
      return NULL;
    }

    bool
    checkSpeciesIsKnown(const std::string& speciesName) const
    {
      return !(theSpeciesListCatalog.find( &speciesName ) == theSpeciesListCatalog.end());
	       
    }

    SpeciesListCatalog& 
    getSpeciesCatalog()
    {
      return theSpeciesListCatalog;
    }
      

    const SpeciesListCatalog& 
    getSpeciesCatalog() const
    {
      return theSpeciesListCatalog;
    }

    const ReactionList&
    getReactionList() const
    {
      return theCompleteReactionList;
    }

    // These two functions are the interface that reaction generators use 
    // to record their species.  They add the (speciesName, speciesPtr)
    // and (reactionName, reactionPtr) entries to the catalog respectively.
    // If the object being recorded is new, we also record it as a hit in the
    // total number of species/reactions as well as the delta number of species/
    // reactions.
      
    bool 
    recordSpecies( SpeciesTypePtr pSpecies)
    {
      std::string* pSpeciesName = new std::string( pSpecies->getName() );

      // This is for debugging.
      std::cout << *pSpeciesName << std::endl;
      std::string* npSpeciesName = new std::string( pSpecies->getName() );
      delete npSpeciesName;


      if (theSpeciesListCatalog.find( pSpeciesName) == theSpeciesListCatalog.end() )
	{
                
	  theSpeciesListCatalog.insert( std::make_pair( pSpeciesName, pSpecies) );
	  theDeltaSpeciesList.push_back( pSpecies );
	  return true;
	}
      else
	{
	  delete pSpeciesName;
	  return false;
	}

    }

      

    // You don't own this.  Don't delete it.  
    ReactionTypePtr
    calculateReactionBetweenSpecies( const speciesT& firstSpecies,
				     const speciesT& secondSpecies)
    {
      // If there is no reaction

    }

    bool
    recordReaction( ReactionTypePtr pRxn )
    {
        // In the current form of the application, this check is not necessary, as any reaction
        // is created at most once.

//         if (theCompleteReactionList.end() != std::find(theCompleteReactionList.begin(),
//                                                        theCompleteReactionList.end(),
//                                                        pRxn))
//         {
//             return false;
            

//         }

        if(!pRxn->isStandardReaction())
        {
            std::cerr<< "A reaction is nonstandard." << std::endl;
            throw 666;
        }

      theCompleteReactionList.push_back( pRxn );
      theDeltaReactionList.push_back( pRxn );

      switch( pRxn->getReactants().size() )
      {
          SpeciesTypePtr theSpeciesPtr;
      case 0:
          // Register in the map of creation reactions
          zeroSubstrateRxns.push_back( pRxn );
          break;
      case 1:
          // Register in a multi map of 1->? reactions
          
          theSpeciesPtr = pRxn->getReactants().begin()->first;
          singleSubstrateRxns.insert( std::make_pair(theSpeciesPtr, pRxn ) );
          break;

      case 2:

          BOOST_FOREACH( typename fnd::basicReaction<SpeciesType>::multMap::value_type vt, 
                         pRxn->getReactants() )
          {
              theSpeciesPtr = vt.first;
              doubleSubstrateRxns.insert( std::make_pair( theSpeciesPtr, pRxn) );
          }

          break;

      default:
          for(unsigned int i = 0; i != 6000; ++i)
          {
              // Spider Jerusalem.
              std::cerr << "FUCK " << std::endl;
          }

          throw 666;
          break;
      }

      return true;
    }

      
    int 
    getPopulation( const std::string& rName) const
    {
      return findSpecies(rName).getPop();
    }
        
    // These are not used by anyone at the moment.  It might be better to use them instead
    // of the plain old recordSpecies/recordReaction functions, but who knows.
    void mustRecordSpecies( SpeciesTypePtr pSpecies ) throw(utl::xcpt)
    {
      if ( !recordSpecies( pSpecies ) )
	{
	  throw DuplicatedCatalogEntryXcpt( pSpecies->getName() );
	}
    }
  
    void
    mustRecordReaction( ReactionTypePtr pRxn) throw(utl::xcpt)
    {
      if ( recordReactions( pRxn ) ) 
	{
	  return;
	}
      else
	{
	  throw DuplicatedCatalogEntryXcpt( pRxn->getName() );
	}

    }

    unsigned int 
    getTotalNumberSpecies() const
    {
      return theSpeciesListCatalog.size();
    }

    unsigned int 
    getTotalNumberReactions() const
    {
      return theCompleteReactionList.size();
    }

    // The getDeltaNumber(Species|Reactions) functions return the number
    // of those objects created since the last time the recordCurrentState
    // funsource ction was source called.
    unsigned int 
    getNumberDeltaReactions() const
    {
      return theDeltaReactionList.size();
    }

    unsigned int 
    getNumberDeltaSpecies() const
    {
      return theDeltaSpeciesList.size();
    }

    void resetCurrentState()
    { 
      theDeltaSpeciesList.clear();
      theDeltaReactionList.clear();
    }

    void incrementNetworkBySpeciesName(const std::string& rName) throw(utl::xcpt)
    {
      SpeciesCatalogCIter iter = theSpeciesListCatalog.find( &rName );

      if (iter != theSpeciesListCatalog.end() )
	{
	  iter->second->expandReactionNetwork();
	}
      else
	{
	  throw NoSuchSpeciesXcpt( rName );
	}
    }

    ~ReactionNetworkDescription()
    {
      // FIXME...
      // We don't memory manage any SpeciesType* or ReactionType*, but we do memory 
      // manage the string* in theSpeciesListCatalog.


      BOOST_FOREACH( typename SpeciesCatalog::value_type i, theSpeciesListCatalog)
	{
	  delete i.first;
	}
    }

    ReactionNetworkDescription(bool reserveMemory = true)
    {
      if (reserveMemory)
	{
	  // theSpeciesListCatalog.reserve( getPredictedSpeciesNetworkSize() );
	  // theCompleteReactionList.reserve( getPredictedSpeciesNetworkSize() );
	  // static_cast<unsigned int>(getPredictedSpeciesReactionRatio() + 1) );
	}
    }

      
      


    // -- The pointers to the strings in theSpeciesListCatalog ARE memory managed.
    // -- The pointers to the species and reactions ARE NOT memory managed here.

    SpeciesCatalog theSpeciesListCatalog;
    // SpeciesListCatalog theSpeciesListCatalog;
    ReactionList theCompleteReactionList;
    ReactionList unaryReactionList;
    ReactionList binaryReactionList;
        
    SpeciesList    theDeltaSpeciesList;
    ReactionList    theDeltaReactionList;

    // 0->1, 1->0, 1->1, 1->2, 2->1

      ReactionList zeroSubstrateRxns;
      ParticipatingSpeciesRxnMap singleSubstrateRxns;
      ParticipatingSpeciesRxnMap doubleSubstrateRxns;
      

  };    


    
}
#endif // RXNNETWORKCATALOG_HH
    
