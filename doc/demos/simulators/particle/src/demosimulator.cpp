#include "demosimulator.hpp"
#include "utl/utility.hh"
#include <fstream>
#include <boost/tokenizer.hpp>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

void SimpleSimulator::createModelFromFile(const std::string& modelFile, std::map<std::string, int>& model)
{
    model.clear();

    std::ifstream inputfile;
    inputfile.open( modelFile.c_str() );

    std::string current_line;

    std::cout << "#######################################################" << "\n" 
              << "Parsing '" << modelFile << "'\n\n" << std::endl;
    while( std::getline(inputfile, current_line))
    {
        std::vector<std::string> tokenVector;

        utl::tokenize( current_line, tokenVector, " ");
        if (tokenVector.size() == 0 ) continue;


        if(tokenVector.size() != 2)
        {
            if (tokenVector.size() == 1 && tokenVector[0] == "") continue;

            std::cerr << "Error, bad line '" << current_line << "'.\n";
            int num = 1;
            BOOST_FOREACH(std::string str, tokenVector)
            {
                std::cerr << num++ << ":\t" << str << '\n';
            }
            assert(tokenVector.size() == 2);
            
        }

        std::string name = tokenVector[0];
        std::string numAsString = tokenVector[1];
        int population;
        
        utl::from_string(population, numAsString);
        std::cout << "Parsing: " << name << " -> " << population << std::endl;

        model.insert( std::make_pair( name, population) );
    }
}

SimpleSimulator::SimpleSimulator( std::string rulesfile,
                                  std::string modelfile )
{
    srand(time(NULL) );
    attachRuleFile( rulesfile);
    attachModelFile( modelfile );
}

void SimpleSimulator::attachRuleFile(std::string rulesfile)
{
    if (speciesReactionGenerator.getModelHasBeenLoaded() )
    {
        std::cerr << "Modelfile has already been loaded.  Error." << std::endl;
        exit(1);
    }
    
    speciesReactionGenerator.attachFileName( rulesfile );
}


void SimpleSimulator::attachModelFile( std::string modelfile)
{

    theModel.clear();

    std::map<std::string, int> theNewModel;
    createModelFromFile( modelfile, theNewModel );

    if (!assertModelValidity( theNewModel ) )
    {
        std::cerr << "Error, the model is not valid.  Crashing." << std::endl;
        exit( 1 );
    }

    theModel.swap( theNewModel );
    engageModel();
}


bool SimpleSimulator::assertModelValidity(const std::map<std::string, int>& model)
{
    if(model.size() == 0)
    {
        cerr << "Error.  No species are present in the model." << endl;
        return false;
    }

    BOOST_FOREACH( const modelPairType& thePair, model)
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

void SimpleSimulator::engageModel()
{
    BOOST_FOREACH( const modelPairType& thePair, theModel)
    {
        if (thePair.second > 0)
        {
            speciesReactionGenerator.incrementNetworkBySpeciesName( thePair.first );
        }
    }
}


void SimpleSimulator::executeReaction( mzr::moleculizer::ReactionTypePtr ptrRxn)
{

    BOOST_FOREACH(const mzr::moleculizer::ReactionType::multMap::value_type& vt, ptrRxn->getReactants())
    {
        std::cout << '-' << vt.second << ' ' << vt.first->getName() << '\n';

        std::string name = vt.first->getName();
        theModel[ name ] -= vt.second;
    }

    BOOST_FOREACH( const mzr::moleculizer::ReactionType::multMap::value_type& vt, ptrRxn->getProducts())
    {
        std::cout << '+' << vt.second << ' ' << vt.first->getName() << '\n';

        std::string name = vt.first->getName();
        theModel[ name ] += vt.second;
    }

    std::cout << std::endl;

}


void SimpleSimulator::printState() const
{
    typedef std::pair<std::string, int> entry;
    BOOST_FOREACH(const entry& ent, theModel)
    {
        std::cout << "\t" << ent.first << ":\t" << ent.second << std::endl;
    }
}
