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
#include "noSuchNameManglerXcpt.hh"
#include "nmrManglerFactory.hh"
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

    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;

  public:
    nmrUnit(mzr::moleculizer& rMoleculizer,
            mzr::mzrUnit& refMzrUnit,
            bnd::molUnit& refMolUnit,
            plx::plexUnit& refPlexUnit)
      :
      unit("nmr", 
           rMoleculizer),
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      ptrNameManglerFactory( new NameManglerFactory<molT> ),
      ptrNameAssembler( NULL )
    {
      setDefaultNameMangler( manglernames::compactManglerName );
    }
    
    ~nmrUnit()
    {
      delete ptrNameManglerFactory;
      delete ptrNameAssembler;
    }

    plx::mzrPlexSpecies* 
    makePlexFromName(const std::string mangledName)
    {
      return NULL;
    }
    
    const NameAssembler<molT>* getNameAssembler() const throw( utl::xcpt );
    void setDefaultNameMangler( const std::string& nameManglerName) throw( NoSuchNameManglerXcpt  );
    

    mzr::mzrSpecies*
    getSpeciesFromName( const std::string& speciesName)
    {
      return NULL;
    }

    virtual void
    parseDomInput(xmlpp::Element* pRootElt,
		  xmlpp::Element* pModelElt,
		  xmlpp::Element* pStreamsElt) throw(std::exception)
    {
      try
        {

          // Redo this, as it is much (maybe) worse than it should be.
          std::string namingConventionStyle = utl::dom::mustGetAttrString( pModelElt,
                                                                           eltName::namingConvention);
          setDefaultNameMangler( namingConventionStyle );

        }
      catch(...)
        {
          return;
        }
    }

    virtual void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
    {
      // For now, do nothing.  It may be best that when the parseDomInput function is 
      // implemented, we will want to insert the generator we are using too.
    }
    
  private:

    NameManglerFactory<molT>* ptrNameManglerFactory;
    NameAssembler<molT>* ptrNameAssembler;
    
  };

}


#include "nmrUnitImpl.hh"

#endif

