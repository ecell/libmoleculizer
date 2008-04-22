/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008  Nathan Addy
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

#include "mzr/mzrSpecies.hh"
#include "nmrUnit.hh"
#include "plex/plexUnit.hh"
#include <string>

namespace nmr
{

    const NameAssembler*
    nmrUnit::getNameEncoder() const 
        throw(utl::xcpt)
    {
        if (!ptrNameAssembler)
        {
            throw utl::xcpt("Error in nmrUnit::getNameEncoder.  No ptrNameAssembler set yet!!!.");
        }

        return ptrNameAssembler;
    }
    
    mzr::mzrSpecies*
    nmrUnit::getSpeciesFromName( const std::string& speciesName)
    {
        if(0)
        {
            // TODO: If the name is present in moleculizer's list, return that.  
        }
        else
        {
            ComplexOutputState aCOS = getNameEncoder()->createOutputStateFromName( speciesName );
            ComplexSpecies theNewComplexSpecies( aCOS );

            mzr::mzrSpecies* newMzrSpecies = pPlexUnit->constructNewPlexSpeciesFromComplexSpecies( theNewComplexSpecies );
            return newMzrSpecies;
        }

        return NULL;
    }

    void
    nmrUnit::setDefaultNameEncoder( const std::string& nameEncoderName) throw( NoSuchNameEncoderXcpt)
    {
        delete ptrNameAssembler;
        ptrNameAssembler = ptrNameEncoderFactory->create(nameEncoderName);
    }



    void nmrUnit::parseDomInput(xmlpp::Element* pRootElt, xmlpp::Element* pModelElt, xmlpp::Element* pStreamsElt) 
        throw(std::exception)
    {
        // TODO: write me
    }

    void nmrUnit::insertStateElts(xmlpp::Element* pRootElt) 
        throw(std::exception)
    {
        // TODO: write me.
    }
}
