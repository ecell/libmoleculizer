/////////////////////////////////////////////////////////////////////////////
// libMoleculizer - a species and reaction generator for reaction networks
// Copyright (C) 2008 The Molecular Sciences Institute
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


#include <iostream>
#include "demostochasticsimulator.hpp"
#include <algorithm>

SimpleStochasticSimulator::SimpleStochasticSimulator(std::string rulesfile, std::string modelfile)
    :
    SimpleSimulator(rulesfile, modelfile)
{

    recordNewReactions();

}

void SimpleStochasticSimulator::recordNewReactions() 
{
    std::copy( speciesReactionGenerator.theDeltaReactionList.begin(),
               speciesReactionGenerator.theDeltaReactionList.end(),
               std::back_inserter( reactions ) );

    std::cout << "Adding: " << std::endl;
    BOOST_FOREACH(mzr::mzrReaction* rxn, speciesReactionGenerator.theDeltaReactionList)
    {
        printRxn( rxn );
    }

    assert( reactions.size() == speciesReactionGenerator.theCompleteReactionList.size() );

    speciesReactionGenerator.printAll();


//     assert( reactions.size() == speciesReactionGenerator.zeroSubstrateRxns.size() + \
//             speciesReactionGenerator.singleSubstrateRxns.size() +    \
//             speciesReactionGenerator.doubleSubstrateRxns.size() );
            
    speciesReactionGenerator.resetCurrentState();
}

mzr::mzrReaction* SimpleStochasticSimulator::calculateReactionToFire()
{

    // Get a list of potential reactions, in this case a list of reactions
    // that have all positive numbers of substrates.

    std::vector<mzr::moleculizer::ReactionTypePtr> potentialReactions;
    getReactionsWithPositivePropensity( potentialReactions );

    if (potentialReactions.size() == 0)
    {
        std::cout << "There are no reactions with positive propensity..." << std::endl;
        return NULL;
    }


    // Pick a random one. 
    int randomIndex = rand() % potentialReactions.size();
    mzr::moleculizer::ReactionTypePtr randomReactionPtr = potentialReactions[ randomIndex ];

    return randomReactionPtr;
}


void SimpleStochasticSimulator::getReactionsWithPositivePropensity(std::vector<mzr::mzrReaction*>& okReactions)
{
    okReactions.clear();

    BOOST_FOREACH(mzr::mzrReaction* rxnptr, reactions)
    {
        if( reactionHasPositiveSubstrates( rxnptr ) )
        {
            okReactions.push_back( rxnptr );
        }
    }
}

bool SimpleStochasticSimulator::reactionHasPositiveSubstrates(const mzr::mzrReaction* rxnPtr)
{
    typedef std::pair<mzr::mzrSpecies*, int> PairType;
    BOOST_FOREACH( const PairType& pr, rxnPtr->getReactants() )
    {
        std::string substrateName( pr.first->getName() );

        if (theModel.find( substrateName) == theModel.end() ) 
        {
            std::cerr << "Error." << std::endl;
        }

        if (theModel[ substrateName ] < pr.second ) return false;
    }

    return true;
    
}

void SimpleStochasticSimulator::singleStep() 
{
    mzr::mzrReaction* reaction = calculateReactionToFire();

    if (reaction)
    {
        executeReaction( reaction );
    }

    recordNewReactions();

    return;
}


