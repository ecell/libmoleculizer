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


#include "ReactionNetworkGenerator.hpp"
#include "mzr/unitsMgr.hh"
#include "nmr/nmrUnit.hh"


Reaction::Reaction(const fnd::basicReaction<mzr::mzrSpecies>& aReaction)
    :
    rate( aReaction.getRate() )
{
    for (CoreRxnType::multMap::const_iterator const_iter = aReaction.getReactants().begin();
         const_iter != aReaction.getReactants().end();
         ++const_iter)
    {
        for(unsigned int ii = 0; ii != const_iter->second; ++ii)
        {
            substrates.push_back( *const_iter->first );
        }
    }

    for (CoreRxnType::multMap::const_iterator const_iter = aReaction.getProducts().begin();
         const_iter != aReaction.getProducts().end();
         ++const_iter)
    {
        for(unsigned int ii = 0; ii != const_iter->second; ++ii)
        {
            products.push_back( *const_iter->first );
        }
        
    }
}

std::vector<Reaction>
ReactionNetworkGenerator::getBinaryReactions(const std::string& species1,
                                             const std::string& species2) throw( mzr::IllegalNameXcpt )
{
  mzr::mzrSpecies* pSpeciesOne = ptrMoleculizer->getSpeciesWithName( species1 );
  mzr::mzrSpecies* pSpeciesTwo = ptrMoleculizer->getSpeciesWithName( species2 );
    
  mzr::mzrReaction* reactionBetweenSpecies = 
    ptrMoleculizer->findReactionWithSubstrates(pSpeciesOne, pSpeciesTwo);

  std::vector<Reaction> theVector;
  
  if( reactionBetweenSpecies )
    {
      theVector.push_back( Reaction(*reactionBetweenSpecies));
    }
  
  return theVector;
}


std::vector<Reaction>
ReactionNetworkGenerator::getUnaryReactions(const std::string& species1) throw( mzr::IllegalNameXcpt )
{
  mzr::mzrSpecies* pSpeciesOne = ptrMoleculizer->getSpeciesWithName( species1);

  std::vector<mzr::mzrReaction*> aVector;
  
  ptrMoleculizer->findReactionWithSubstrates(pSpeciesOne, aVector);
  
  std::vector<Reaction> theVector;
  theVector.reserve( aVector.size() );

  
  BOOST_FOREACH( mzr::mzrReaction* ptrMzrRxn, aVector)
    {
      theVector.push_back( Reaction(*ptrMzrRxn));
    }
  
  return theVector;
}

Species 
ReactionNetworkGenerator::getSpecies(const std::string& speciesName) throw( mzr::IllegalNameXcpt )
{
    return Species( *ptrMoleculizer->getSpeciesWithName( speciesName ) );
}


bool
ReactionNetworkGenerator::checkSpeciesNameLegality(const std::string& speciesOne) throw (mzr::IllegalNameXcpt )
{
    try
    {
        ptrMoleculizer->getSpeciesWithName( speciesOne );
        return true;
    }
    catch(...)
    {
      return false;
    }
}
