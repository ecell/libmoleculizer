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

#ifndef REACTIONNETWORKGENERATOR_HPP
#define REACTIONNETWORKGENERATOR_HPP

#include "mzr/moleculizer.hh"
#include "mzr/mzrException.hh"
#include "fnd/basicReaction.hh"
#include "mzr/mzrSpecies.hh"
#include <vector>
#include <string>

class Species
{
public:
    Species()
        :
        name(""),
        mass(0.0f){}

    Species(const mzr::mzrSpecies& aMzrSpecies)
        :
        name( aMzrSpecies.getName() ),
        mass( aMzrSpecies.getWeight() )
    {}

    std::string getName() const { return name;}
    float getMass() const {return mass;}
    
private:
    std::string name;
    float mass;
};

class Reaction
{
public:
    typedef fnd::basicReaction<mzr::mzrSpecies> CoreRxnType;

    Reaction(){}
    Reaction(const CoreRxnType& aReaction);

    float getRate() const
    {
        return rate;
    }

    std::vector<Species> 
    getSubstrates() const
    {
        return substrates;
    }

    std::vector<Species> 
    getProducts() const
    {
        return products;
    }

protected:
    float rate;

    std::vector<Species> substrates;
    std::vector<Species> products;
};


class ReactionNetworkGenerator
{

public:
    ReactionNetworkGenerator()
    {
        ptrMoleculizer = new mzr::moleculizer;
    }

    ~ReactionNetworkGenerator()
    {
        delete ptrMoleculizer;
    }

    void addRules(const std::string& filename) throw( mzr::BadRulesDefinitionXcpt )
    {
        try
        {
            delete ptrMoleculizer;
            ptrMoleculizer = new mzr::moleculizer;
            ptrMoleculizer->attachFileName( filename );
        }
        catch(mzr::BadRulesDefinitionXcpt e)
        {
            e.warn();

            if(ptrMoleculizer)
            {
                delete ptrMoleculizer;
                ptrMoleculizer = NULL;
            }

            throw e;
        }
        catch(utl::xcpt x)
        {
            delete ptrMoleculizer;
            x.warn();
            throw mzr::BadRulesDefinitionXcpt();
        }
        catch(...)
        {
            delete ptrMoleculizer;
            utl::xcpt x("Unknown Error while adding rules. (FDLLAE)");
            x.wailAndBail();
        }
    }
    
    std::vector<Reaction>
    getBinaryReactions(const std::string& species1,
                       const std::string& species2) throw( mzr::IllegalNameXcpt );

    std::vector<Reaction>
    getUnaryReactions(const std::string& species1) throw (mzr::IllegalNameXcpt );

    Species 
    getSpecies(const std::string& species)throw (mzr::IllegalNameXcpt );

    bool 
    checkSpeciesNameLegality( const std::string& species1) throw (mzr::IllegalNameXcpt );

    void 
    showAllSpecies() const
    {
        ptrMoleculizer->DEBUG_showAllSpecies();
    }


protected:
    mzr::moleculizer* ptrMoleculizer;
};


#endif
