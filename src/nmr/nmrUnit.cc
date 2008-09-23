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
// along with Moleculizer; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
/////////////////////////////////////////////////////////////////////////////

#include "plex/mzrPlex.hh"
#include "mzr/mzrSpecies.hh"
#include "plex/mzrPlexSpecies.hh"
#include "nmrUnit.hh"
#include "nmr/nmrEltName.hh"
#include "plex/plexUnit.hh"
#include "nmrExceptions.hh"
#include "mol/mzrMol.hh"

#include <string>

namespace nmr
{
    const NameAssembler*
    nmrUnit::getNameEncoder() const 
        throw( MissingNameEncoderXcpt )
    {
        if (!ptrNameAssembler)
        {
            throw MissingNameEncoderXcpt();
        }

        return ptrNameAssembler;
    }
    
    mzr::mzrSpecies*
    nmrUnit::constructSpeciesFromName( const std::string& encodedSpeciesName) 
    {
      try
	{
	  // 1.  The currently set nameEncoder has the responsibility of decoding the thing into a complexOutputState.
	  ComplexOutputState aCOS = getNameEncoder()->createOutputStateFromName( encodedSpeciesName );

	  // TODO: Ensure that this does  not register the constructed plexSpecies into 
	  plx::mzrPlexSpecies* newMzrSpecies = pPlexUnit->constructNewPlexSpeciesFromComplexOutputState( aCOS );
	  
	  return newMzrSpecies;
	}
      // I should separate the reasons things went bad here.
      catch(...)
	{
	  throw nmr::IllegalNameXcpt(encodedSpeciesName);
	}

    }

    void
    nmrUnit::setDefaultNameEncoder( const std::string& nameEncoderName) throw( NoSuchNameEncoderXcpt)
    {
        delete ptrNameAssembler;
        ptrNameAssembler = ptrNameEncoderFactory->create(nameEncoderName);
    }



    void 
    nmrUnit::parseDomInput(xmlpp::Element* pRootElt, 
                           xmlpp::Element* pModelElt,
                           xmlpp::Element* pStreamElt) throw(std::exception)
    {
        // TODO: Debug nmrUnit::parseDomInput and make sure it all works.
        // For the moment, the one thing (the naming strategy that should be used) is contained 
        // in an attribute entitled "naming-convention".  The value of that string should 
        // be passed to nmrUnit::setDefaultNameEncoder.  


        try
        {
            std::string strategy = utl::dom::mustGetAttrString(pModelElt, eltName::namingConvention);
            setDefaultNameEncoder(strategy);
        }
        catch(NoSuchNameEncoderXcpt xcpt)
        {
            xcpt.wailAndBail();
        }
        catch(utl::xcpt x)
        {
            // This is somewhat bad mojo, but it would be if something is an utl::xcpt
            // but has not already been caught.
            
            // Do nothing because this element is not mandatory.
            x.warn();
        }
        
    }

    void nmrUnit::insertStateElts(xmlpp::Element* pRootElt) 
        throw(std::exception)
    {
        // TODO: Debug nmrUnit::insertStateElts and make sure it all works.
        xmlpp::Element* modelElt = utl::dom::mustGetUniqueChild(pRootElt, 
                                                                mzr::eltName::model);
        
        modelElt->set_attribute( eltName::namingConvention,
                                getNameEncoder()->getName() );
    }
}
