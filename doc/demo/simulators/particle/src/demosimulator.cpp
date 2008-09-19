#include "demosimulator.hpp"
#include "utl/utility.hh"
#include <fstream>
#include <boost/tokenizer.hpp>

void SimplifiedStochasticTypeSimulator::createModelFromFile(const std::string& modelFile, std::map<std::string, int>& model)
{
    model.clear();

    std::ifstream inputfile;
    inputfile.open( modelFile.c_str() );

    std::string current_line;
    while( std::getline(inputfile, current_line))
    {
        boost::tokenizer<> tok(current_line);
        std::vector<std::string> tokenVector(tok.begin(), tok.end());
        
        assert(tokenVector.size() == 2);

        std::string name = tokenVector[0];
        std::string numAsString = tokenVector[1];
        int population;
        
        utl::from_string(population, numAsString);

        model.insert( std::make_pair( name, population) );
    }
}

SimplifiedStochasticTypeSimulator::SimplifiedStochasticTypeSimulator(std::string modelfile, 
                                                                     std::string rulesfile)
    :
    speciesReactionGenerator(),
    theModel(),
    theReactionList()
{
    speciesReactionGenerator.attachFileName( rulesfile );
    createModelFromFile( modelfile, theModel);

    if (!assertModelValidity()) 
    {
        std::cerr << "Model is invalid." << std::endl;
    }

    engageModel();
}
bool SimplifiedStochasticTypeSimulator::assertModelValidity()
{
    if(theModel.size() == 0)
    {
        cerr << "Error.  No species are present in the model." << endl;
        return false;
    }

    BOOST_FOREACH( const modelPairType& thePair, theModel)
    {
        // Find out if this is a legal name.
        try
        {
            speciesReactionGenerator.getSpeciesWithName( thePair.first );
        }
        catch(mzr::IllegalNameXcpt)
        {
            std::cerr << "Illegal name: " << thePair.first << std::endl;
            return false;
        }
        catch(...)
        {
            return false;
        }
    }
    
    return true;
}

void SimplifiedStochasticTypeSimulator::engageModel()
{
    BOOST_FOREACH( const modelPairType& thePair, theModel)
    {
        speciesReactionGenerator.incrementNetworkBySpeciesName( thePair.first );
    }
}

void SimplifiedStochasticTypeSimulator::singleStep()
{
    
}

void SimplifiedStochasticTypeSimulator::attachRuleFile(std::string rulesfile)
{
    speciesReactionGenerator.attachFileName( rulesfile );
}


