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

#include "demosimulator.hpp"
#include  "utl/defs.hh"
#include "utl/utility.hh"
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
using std::cout;
using std::cerr;
using std::endl;


void SimpleSimulator::createModelFromFile (const std::string& modelFile, std::map<std::string, int>& model)
{
    model.clear();

    std::ifstream inputfile;
    inputfile.open ( modelFile.c_str() );

    std::string current_line;
    while ( std::getline (inputfile, current_line) )
    {
        std::vector<std::string> tokenVector;

        utl::tokenize ( current_line, tokenVector, " ");
        if (tokenVector.size() == 0 ) continue;


        if (tokenVector.size() != 2)
        {
            if (tokenVector.size() == 1 && tokenVector[0] == "") continue;

            std::cerr << "Error, bad line '" << current_line << "'.\n";
            int num = 1;
            BOOST_FOREACH (std::string str, tokenVector)
            {
                std::cerr << num++ << ":\t" << str << '\n';
            }
            assert (tokenVector.size() == 2);

        }

        std::string name = tokenVector[0];
        std::string numAsString = tokenVector[1];
        int population;

        utl::from_string (population, numAsString);
        model.insert ( std::make_pair ( name, population) );
    }
}

SimpleSimulator::SimpleSimulator ( std::string rulesfile,
                                   std::string modelfile )
{
    srand (time (NULL) );
    attachRuleFile ( rulesfile);
    
    std::cout << "SS: State after attaching rules, before attaching Model:" << std::endl;
    speciesReactionGenerator.printAll();

    std::cout << "---------(END)" << std::endl;
    
    attachModelFile ( modelfile );

    std::cout << "---------(BEGIN)" << std::endl;

    std::cout << "SS: State after attaching rules, model" << std::endl;
    speciesReactionGenerator.printAll();

    std::cout << "---------(END)" << std::endl;
}

void SimpleSimulator::attachRuleFile (std::string rulesfile)
{
    if (speciesReactionGenerator.getModelHasBeenLoaded() )
    {
        std::cerr << "Modelfile has already been loaded.  Error." << std::endl;
        exit (1);
    }

    speciesReactionGenerator.attachFileName ( rulesfile );
}


void SimpleSimulator::attachModelFile ( std::string modelfile)
{

    theModel.clear();

    std::map<std::string, int> theNewModel;
    createModelFromFile ( modelfile, theNewModel );

    if (!assertModelValidity ( theNewModel ) )
    {
        std::cerr << "Error, the model is not valid.  Crashing." << std::endl;
        exit ( 1 );
    }

    theModel.swap ( theNewModel );
    engageModel();
}


bool SimpleSimulator::assertModelValidity (const std::map<std::string, int>& model)
{
    if (model.size() == 0)
    {
        cerr << "Error.  No species are present in the model." << endl;
        return false;
    }

    BOOST_FOREACH ( const modelPairType& thePair, model)
    {
        // Find out if this is a legal name.
        try
        {
            speciesReactionGenerator.getSpeciesWithName ( thePair.first );
        }
        catch (mzr::IllegalNameXcpt)
        {
            std::cerr << "Illegal name: " << thePair.first << std::endl;
            return false;
        }
        catch (...)
        {
            return false;
        }
    }

    return true;
}

void SimpleSimulator::engageModel()
{
    BOOST_FOREACH ( const modelPairType& thePair, theModel)
    {
        if (thePair.second > 0)
        {
            speciesReactionGenerator.incrementNetworkBySpeciesName ( thePair.first );
        }
    }
}


void SimpleSimulator::executeReaction ( mzr::moleculizer::ReactionTypePtr ptrRxn)
{
    std::cout << "Executing: " << ptrRxn->getName() << std::endl;

    BOOST_FOREACH (const mzr::moleculizer::ReactionType::multMap::value_type& vt, ptrRxn->getReactants() )
    {
        std::cout << '-' << vt.second << ' ' << vt.first->getName() << '\n';

        std::string name = vt.first->getName();
        theModel[ name ] -= vt.second;
    }

    BOOST_FOREACH ( const mzr::moleculizer::ReactionType::multMap::value_type& vt, ptrRxn->getProducts() )
    {
        std::cout << '+' << vt.second << ' ' << vt.first->getName() << '\n';

        std::string name = vt.first->getName();

        if (theModel.find ( name ) == theModel.end() )
        {
            std::cout << "Creating new species '" << name << "'.";
            theModel[ name ] = vt.second;
        }
        else
        {
            theModel[ name ] += vt.second;
        }
    }

    std::cout << std::endl;

}


void SimpleSimulator::printState() const
{
    std::cout << "SimpleSimulator state:" << std::endl;
    typedef std::pair<std::string, int> entry;
    BOOST_FOREACH (const entry& ent, theModel)
    {
        std::cout << "\t" << ent.first << ":\t" << ent.second << std::endl;
    }
    std::cout << "Moleculizer State" << std::endl;

    speciesReactionGenerator.printAll();
}
