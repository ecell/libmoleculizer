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


#include "demoparticlesimulator.hpp"

void SimpleParticleSimulator::singleStep()
{
        doSingleUnaryReaction();

    std::vector<std::string> stringPtrVector;

    BOOST_FOREACH( const modelPairType& mpt, theModel)
    {
        if (mpt.second > 0)
        {
            stringPtrVector.push_back( mpt.first );
        }
    }

    // Pick two particles at random.

    std::string speciesNameOne, speciesNameTwo;
    bool illegalpick = true;
    
    do
    {
        int firstIndex = rand() % stringPtrVector.size();
        int secondIndex = rand() % stringPtrVector.size();
        speciesNameOne = stringPtrVector[ firstIndex ];
        speciesNameTwo = stringPtrVector[ secondIndex ];

        if (speciesNameOne == speciesNameTwo && theModel[speciesNameOne] >= 2)
        {
            illegalpick = false;
        }
        else
        {
            illegalpick = false;
        }
    }
    while(illegalpick);


    std::cout << "Collision between: " << '\n'
              << '\t' << speciesNameOne << '\n'
              << '\t' << speciesNameTwo << std::endl;

    // Look up these species in the moleculizer list.  
    mzr::moleculizer::SpeciesTypePtr speciesPtrOne = speciesReactionGenerator.getSpeciesWithName( speciesNameOne );
    mzr::moleculizer::SpeciesTypePtr speciesPtrTwo = speciesReactionGenerator.getSpeciesWithName( speciesNameTwo );

    std::vector<mzr::moleculizer::ReactionTypePtr> reactionVector;
    speciesReactionGenerator.findReactionWithSubstrates( speciesPtrOne,
                                                         speciesPtrTwo,
                                                         reactionVector );

    if (reactionVector.size() == 0)
    {
        std::cout << "\tNORES \t-- (No reactions between)." << std::endl << std::endl;
        return;
    }

    
    std::cout << "\t!!! reaction." << std::endl;
    executeReaction( reactionVector[ rand() % reactionVector.size() ] );
}


void SimpleParticleSimulator::doSingleUnaryReaction()
{

    std::vector<std::string> stringPtrVector;
    BOOST_FOREACH( const modelPairType& mpt, theModel)
    {
        if (mpt.second > 0)
        {
            stringPtrVector.push_back( mpt.first );
        }
    }

    std::string particleName = stringPtrVector[rand() % stringPtrVector.size() ];

    std::cout << "Possible unary reaction: " << '\n'
              << '\t' << particleName << endl;

    mzr::moleculizer::SpeciesTypePtr particlePtr = speciesReactionGenerator.getSpeciesWithName( particleName );

    std::vector<mzr::moleculizer::ReactionTypePtr> reactionVector;
    speciesReactionGenerator.findReactionWithSubstrates(particlePtr,
                                                        reactionVector);

    if ( reactionVector.size() == 0 )
    {
        std::cout << "\tNORES \t-- (No unary reactions).\n" << std::endl;
        return;
    }

    std::cout << "\t!!! reaction." << std::endl;

    mzr::moleculizer::ReactionTypePtr rxn = reactionVector[ rand() % reactionVector.size() ];

    executeReaction( rxn );

}
