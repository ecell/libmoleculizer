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
