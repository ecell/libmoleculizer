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

#ifndef NMRUNIT_HH
#define NMRUNIT_HH

#include "nmrEltName.hh"
#include "nameAssembler.hh"
#include "nameEncoderFactory.hh"
#include "mzr/mzrSpecies.hh"
#include "plex/mzrPlexSpecies.hh"
#include "mzr/unit.hh"
#include "mol/molUnit.hh"
#include "plex/plexUnit.hh"

namespace nmr
{

  template <typename molT>
  class nmrUnit : public mzr::unit
  {
    typedef molT MolType;

  public:
    nmrUnit(mzr::moleculizer& rMoleculizer)
      :
      unit("nmr", 
           rMoleculizer),
      pMzrUnit(NULL),
      pMolUnit(NULL),
      pPlexUnit(NULL),
      ptrNameEncoderFactory( new NameEncoderFactory<molT> ),
      ptrNameAssembler( NULL )
    {
      setDefaultNameEncoder( manglernames::compactEncoderName );
    }
    
    ~nmrUnit()
    {
      // Don't delete any pointers to Units.
      delete ptrNameEncoderFactory;
      delete ptrNameAssembler;
    }

    plx::mzrPlexSpecies* 
    makePlexFromName(const std::string& mangledName) const
    {
      return NULL;
    }
    


    mzr::mzrSpecies*
    getSpeciesFromName( const std::string& speciesName)
    {
      return NULL;
    }
    
    const NameAssembler<molT>* 
    getNameEncoder() const throw( utl::xcpt );

    void 
    setDefaultNameEncoder( const std::string& nameEncoderName) throw( NoSuchNameEncoderXcpt  );

    void
    setMzrUnit(mzr::mzrUnit* ptrMzrUnit)
    {
      pMzrUnit = ptrMzrUnit;
    }

    void 
    setPlexUnit(plx::plexUnit* ptrPlexUnit)
    {
      pPlexUnit = ptrPlexUnit;
    }

    void setMolUnit(bnd::molUnit* ptrMolUnit)
    {
      pMolUnit = ptrMolUnit;
    }

    void parseDomInput(xmlpp::Element* pRootElt, xmlpp::Element* pModelElt, xmlpp::Element* pStreamsElt) throw(std::exception);
    void insertStateElts(xmlpp::Element* pRootElt) throw(std::exception);

  private:

    mzr::mzrUnit* pMzrUnit;
    bnd::molUnit* pMolUnit;
    plx::plexUnit* pPlexUnit;

    NameEncoderFactory<molT>* ptrNameEncoderFactory;
    NameAssembler<molT>* ptrNameAssembler;
  };
}


#include "nmrUnitImpl.hh"

#endif

