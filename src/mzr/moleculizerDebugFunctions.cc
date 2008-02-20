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
#include "mzr/debug.hh"
#include "utl/string.hh"
#include <algorithm>

namespace mzr
{
    class printObject
    {
    public:
        template <typename objType>
        void operator()(objType* ptrObject)
        {
            cout << ptrObject->getName()  << endl;
        }
    };


    void 
    moleculizer::RunInteractiveDebugMode()
    {
        unsigned int result;
        bool cont = true;

        while( cont )
        {
            std::cout << "\n\n";
            std::cout << "0: Quit" << std::endl;
            std::cout << "1: showNumberSpecies()" << std::endl;
            std::cout << "2: showNumberReactions()" << std::endl;
            std::cout << "3: showSpecies()" << std::endl;
            std::cout << "4: showReactions()" << std::endl;
            std::cout << "5: showNewlyCreated()" << std::endl;
            std::cout << "6: showDeltaSpecies()" << std::endl;
            std::cout << "7: showDeltaReactions()" << std::endl;
            std::cout << "8: showLiveSpecies()" << std::endl;
            std::cout << "9: incrementSpecies()" << endl;
            std::cout << "Give input:\t" << std::endl;
            std::cin >> result;
            std::cout << "\n" << std::endl;
            switch (result)
            {
            case 0:
                cont = false;
                break;
            case 1:
                DEBUG_showNumberSpecies();
                break;
            case 2:
                DEBUG_showNumberReactions();
                break;
            case 3:
                DEBUG_showSpecies();
                break;
            case 4:
                DEBUG_showReactions();
                break;
            case 5:
                DEBUG_showNewlyCreated();
                break;
            case 6:
                DEBUG_showDeltaSpecies();
                break;
            case 7:
                DEBUG_showDeltaReactions();
                break;
            case 8:
                DEBUG_showLiveSpecies();
                break;
            case 9:
                DEBUG_incrementSpecies();
                break;
                    
                
            default:
                std::cout << "Option " << result << " not yet implemented" << std::endl;
                continue;
            }
        }
    }

    void moleculizer::RunProfileMode(unsigned int Num_Iterations)
    {
        for( unsigned int i = 0;
             i != Num_Iterations;
             ++i)
        {
	  if (i % (10) == 0 )
            {
	      cout << i * 100.0 / static_cast<float>(Num_Iterations) << "% done..." << endl;
            }

            recordCurrentState();
            std::list<mzrSpecies*>::reverse_iterator riter = find_if( listOfAllSpecies.rbegin(),
                                                                      listOfAllSpecies.rend(),
                                                                      std::not1( std::mem_fun( &mzrSpecies::hasNotified ) ) );

            if (riter == listOfAllSpecies.rend() )
            {
                cout << "Unexpectedly could not find anything... Quitting..." << endl;
                break;
            }
            cout << "Expanding " << (*riter)->getName() << "..." << endl;
            (*riter)->expandReactionNetwork();
            
            cout << "\t" << getDeltaNumberSpecies() << " species generated\n";
            cout << "\t" << getDeltaNumberReactions() << " reactions generated." << endl;

	    DEBUG_showNumberSpecies();
	    DEBUG_showNumberReactions();
	    
	    cout << "######################################################" << endl;

        }

	DEBUG_showSpecies();
	DEBUG_showReactions();
             
               
    }




    void
    moleculizer::DEBUG_showNumberSpecies() const
    {
        cout << "There are " << getTotalNumberSpecies() << " species in the list." << endl;
    }

    void
    moleculizer::DEBUG_showNumberReactions() const
    {
        cout << "There are " << getTotalNumberReactions() << " reactions in the list." << endl;
    }

    void
    moleculizer::DEBUG_showSpecies() const
    {
        unsigned int index = 0;
        for( utl::catalog<mzrSpecies>::const_iterator iter = canonicalCatalogOfSpecies.begin();
	     iter != canonicalCatalogOfSpecies.end();
	     ++iter, ++index)
	  {
	    cout << index << ":\t" << iter->first << endl;
	  }

// 	list<mzrSpecies*>::const_iterator iter = listOfAllSpecies.begin();
//              iter != listOfAllSpecies.end();
//              ++iter, ++i)
//         {
//             cout << i << ":\t" << (*iter)->getName() << endl;
//         }
    }

    void
    moleculizer::DEBUG_showReactions() const
    {
        unsigned int ndx = 0;
        for( utl::catalog<mzrReaction>::const_iterator iter = canonicalCatalogOfRxns.begin();
	     iter != canonicalCatalogOfRxns.end();
	     ++iter, ++ndx)
        {
            cout << ndx << ":\t" << iter->first << endl;
        }
    }

    
    void
    moleculizer::DEBUG_showNewlyCreated() const
    {
        cout << "The new Species: " << endl;
        DEBUG_showDeltaSpecies();
        cout << "The new Reactions: " << endl;
        DEBUG_showDeltaReactions();
    }

    void
    moleculizer::DEBUG_showDeltaSpecies() const
    {
        for( std::list<mzrSpecies*>::const_iterator iter = getLastSpeciesCreationDelta();
             iter != listOfAllSpecies.end();
             ++iter)
        {
            cout << (*iter)->getName() << endl;
        }
        cout << "\tThere are " << getDeltaNumberSpecies() << " delta species." << endl;
    }

    void moleculizer::DEBUG_showDeltaReactions() const
    {

        for(std::list<mzrReaction*>::const_iterator iter = getLastReactionCreationDelta();
            iter != listOfAllReactions.end();
            ++iter)
        {
            cout << (*iter)->getName() << endl;
        }
        cout << "\tThere are " << getDeltaNumberReactions() << " delta reactions." << endl;
    }

    void moleculizer::DEBUG_showLiveSpecies() const
    {
        std::vector<const mzrSpecies*> speciesVector;
        speciesVector.reserve( listOfAllSpecies.size() );

        // Copy any element which does not match the predicate 
        // "the value of that thing ->hadNotified is true".
        std::remove_copy_if(listOfAllSpecies.begin(),
                            listOfAllSpecies.end(),
                            std::back_inserter( speciesVector ),
                            std::mem_fun( &mzrSpecies::hasNotified ) );

        cout << "There are " << speciesVector.size() << " 'live' species." << endl;
        std::for_each( speciesVector.begin(),
                       speciesVector.end(),
                       printObject() );
        
        
    }

    void moleculizer::DEBUG_incrementSpecies()
    {

        recordCurrentState();
        std::string nameToIncrement;
        std::cin >> nameToIncrement;
        
        try
        {
            incrementNetworkBySpeciesName( nameToIncrement );
        }
        catch( utl::xcpt x)
        {
            cout << x.getMessage() << endl;
            cout << "continuing...." << endl;
        }
        
    }

}
