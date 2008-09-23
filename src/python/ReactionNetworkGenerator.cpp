/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008 The Molecular Sciences Institute
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
/////////////////////////////////////////////////////////////////////////////


#include "ReactionNetworkGenerator.hpp"
#include "mzr/unitsMgr.hh"
#include "nmr/nmrUnit.hh"
#include "mol/mzrModMol.hh"
#include <boost/foreach.hpp>


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


void BasicComplexRepresentation::addMolNameToComplex(const std::string& molName)
{
    mols.push_back(molName);
}

void BasicComplexRepresentation::addBindingToComplex(int molIndex,
                                                     const std::string& bindingName1,
                                                     int secondMolNdx,
                                                     const std::string& bindingName2)
{
    bindings.push_back(std::make_pair( std::make_pair(molIndex, bindingName1),
                                       std::make_pair(secondMolNdx, bindingName2)));
}


void BasicComplexRepresentation::addModificationToComplex(int molIndex, 
                                                          const std::string& modificationSiteName,
                                                          const std::string& modificationValue)
{
    modifications.push_back( std::make_pair(molIndex,
                                            std::make_pair(modificationSiteName,
                                                           modificationValue)));
}

std::vector<Reaction>
ReactionNetworkGenerator::getBinaryReactions(const std::string& species1,
                                             const std::string& species2) throw( mzr::IllegalNameXcpt )
{
  mzr::mzrSpecies* pSpeciesOne = ptrMoleculizer->getSpeciesWithName( species1 );
  mzr::mzrSpecies* pSpeciesTwo = ptrMoleculizer->getSpeciesWithName( species2 );

  pSpeciesOne->expandReactionNetwork();
  pSpeciesTwo->expandReactionNetwork();
    
  std::vector<mzr::mzrReaction*> theVector;
  ptrMoleculizer->findReactionWithSubstrates(pSpeciesOne, pSpeciesTwo, theVector);

  std::vector<Reaction> rxnVector;
  BOOST_FOREACH( mzr::mzrReaction* ptrRxn, theVector)
  {
      rxnVector.push_back( Reaction(*ptrRxn));
  }


  return rxnVector;
}

std::vector<Reaction>
ReactionNetworkGenerator::getUnaryReactions(const std::string& species1) throw( mzr::IllegalNameXcpt )
{

    ptrMoleculizer->incrementNetworkBySpeciesName( species1 );
    mzr::mzrSpecies* pSpeciesOne = (*ptrMoleculizer).theSpeciesListCatalog[&species1];
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

void
ReactionNetworkGenerator::showAllReactions()
{
    ptrMoleculizer->DEBUG_showAllReactions();
}

void
ReactionNetworkGenerator::incrementSpecies(const std::string& species)
{
    ptrMoleculizer->incrementNetworkBySpeciesName( species );
}


int 
ReactionNetworkGenerator::getNumberOfSpecies()
{
    return ptrMoleculizer->getTotalNumberSpecies();
}

int 
ReactionNetworkGenerator::getNumberOfReactions()
{
    return ptrMoleculizer->getTotalNumberReactions();
}



std::string 
ReactionNetworkGenerator::generateNameFromBasicComplexRepresentationStrict(const BasicComplexRepresentation& aBCR)
{
    nmr::ComplexSpecies aComplexSpecies;


    for(unsigned int molNdx = 0; molNdx != aBCR.mols.size(); ++molNdx)
    {

        std::string molName = aBCR.mols[molNdx];
        nmr::MinimalMolSharedPtr ptrMol(new nmr::MinimalMol( molName ) );
        aComplexSpecies.addMolToComplex( ptrMol, utl::stringify(molNdx) );

        bnd::mzrMol* ptrMzrMol = ptrMoleculizer->pUserUnits->pMolUnit->mustFindMol( molName );
        for(bnd::mzrMol::bindingSiteIterator iter = ptrMzrMol->getBindingSitesBegin();
            iter != ptrMzrMol->getBindingSitesEnd();
            ++iter)
        {
            std::string bindingSiteName( iter->getName() );
            ptrMol->addNewBindingSite( bindingSiteName );
        }


        const bnd::mzrModMol* ptrMzrModMol = 
            dynamic_cast<const bnd::mzrModMol* >(ptrMzrMol);
        
        if(ptrMzrModMol)
        {
            typedef std::pair<std::string, int> StrIntPair;
            BOOST_FOREACH(const std::string& strRef, ptrMzrModMol->modSiteNames)
            {
                ptrMol->addNewModificationSite( strRef, ptrMzrModMol->getDefaultModNameForSite(strRef) );
            }
        }

        // This isn't the best way to do this, but because nmr::complexSpecies 
        std::vector<BasicComplexRepresentation::ModType> relevantModifications(aBCR.modifications.begin(),
                                                                               aBCR.modifications.end());
        std::vector<BasicComplexRepresentation::ModType>::iterator newEnd = std::remove_if(relevantModifications.begin(),
                                                                                           relevantModifications.end(),
                                                                                           modificationNotOfNdx(molNdx));

        for(std::vector<BasicComplexRepresentation::ModType>::iterator iter = relevantModifications.begin();
            iter != newEnd;
            ++iter)
        {
            ptrMol->updateModificationState(iter->second.first,
                                            iter->second.second);
        }

    }


    BOOST_FOREACH( const BasicComplexRepresentation::BindingType& bt, aBCR.bindings)
    {
        aComplexSpecies.addBindingToComplex( utl::stringify(bt.first.first),
                                             bt.first.second,
                                             utl::stringify(bt.second.first),
                                             bt.second.second);
    }

    const nmr::NameAssembler* pAppNameAssembler( ptrMoleculizer->pUserUnits->pNmrUnit->getNameEncoder() );
    string aComplexSpeciesName( pAppNameAssembler->createCanonicalName(aComplexSpecies) );
    return aComplexSpeciesName;

}


// This 
std::string 
ReactionNetworkGenerator::generateNameFromBasicComplexRepresentation(const BasicComplexRepresentation& aBCR)
{
    nmr::ComplexSpecies aComplexSpecies;

    unsigned int i = 0;
    BOOST_FOREACH(const std::string& molName, aBCR.mols)
    {
        nmr::MinimalMolSharedPtr ptrMol(new nmr::MinimalMol( molName ) );

        aComplexSpecies.addMolToComplex( ptrMol, utl::stringify(i++) );

        bnd::mzrMol* ptrMzrMol = ptrMoleculizer->pUserUnits->pMolUnit->mustFindMol( molName );

        for(bnd::mzrMol::bindingSiteIterator iter = ptrMzrMol->getBindingSitesBegin();
            iter != ptrMzrMol->getBindingSitesEnd();
            ++iter)
        {
            ptrMol->addNewBindingSite( iter->getName() );
        }

// Do this first thing tomorrow...



        // const cpx::modMol<typename plexFamilyT::molType>* aModMol = 
//             dynamic_cast<const cpx::modMol<typename plexFamilyT::molType>* >(pMol);

//         if(aModMol)
//         {
//             // Get the externalized state....
//             const cpx::modMolState& nuMolParam = aModMol->externState( molParams[molNdx] );

//             for(unsigned int ndx = 0;
//                 ndx != aModMol->modSiteNames.size();
//                 ++ndx)
//             {
                
//                 aMol->addNewModificationSite( aModMol->modSiteNames[ndx],
//                                               nuMolParam[ndx]->getName() );
//             }

//         }



    }

    BOOST_FOREACH( const BasicComplexRepresentation::BindingType& bt, aBCR.bindings)
    {
        aComplexSpecies.addBindingToComplex( utl::stringify(bt.first.first),
                                             bt.first.second,
                                             utl::stringify(bt.second.first),
                                             bt.second.second);
    }

    const nmr::NameAssembler* pAppNameAssembler( ptrMoleculizer->pUserUnits->pNmrUnit->getNameEncoder() );
    string aComplexSpeciesName( pAppNameAssembler->createCanonicalName(aComplexSpecies) );
    return aComplexSpeciesName;

}




