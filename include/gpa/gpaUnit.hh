/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001  Walter Lawrence (Larry) Lok.
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
// Contact information:
//   Larry Lok, Research Fellow          Voice: 510-981-8740
//   The Molecular Sciences Institute      Fax: 510-647-0699
//   2168 Shattuck Ave.                  Email: lok@molsci.org
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#ifndef GPAUNIT_H
#define GPAUNIT_H

/*! \defgroup gpaGroup The gpa unit.
  \ingroup unitsGroup
  \brief Provides means to define G-proteins and examine their state. */

/*! \file gpaUnit.hh
  \ingroup gpaGroup
  \brief Adds gpa mols and commands for their allosteric properties. */

#include "mzr/unit.hh"
#include "mzr/mzrEltName.hh"
#include "plex/plexUnit.hh"
#include "plex/plexEltName.hh"
#include "gpa/gpaEltName.hh"

namespace gpa
{
  /*! \ingroup gpaGroup
    \brief Provides mol property of binding to GTP/GDP as does Gpa1.

    This unit provides the special characteristic of G-proteins and enables
    G-proteins to display sites in an allosteric way, depending on whether they
    are bound to GTP or GDP.

    It also provides commands for selecting complexes (or instances of
    subcomplexes) based on the state of G-proteins that they contain.  */
  class gpaUnit : public mzr::unit
  {
  public:
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    stoch::stochUnit& rStochUnit;
    
    gpaUnit(mzr::moleculizer& rMoleculizer,
	    mzr::mzrUnit& refMzrUnit,
	    bnd::molUnit& refMolUnit,
	    plx::plexUnit& refPlexUnit,
	    stoch::stochUnit& refStochUnit) :
      mzr::unit("gpa",
		rMoleculizer),
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      rStochUnit(refStochUnit)
    {
      inputCap.addReactionGenName(eltName::gpaExchangeGen);
      inputCap.addReactionGenName(eltName::gpaRevertGen);

      // Register omniplex nodes for processing by the plexUnit.
      // These are "enabling complexes" for reaction generators.
      const std::string slash("/");
      std::ostringstream gpaExGensXpath;
      gpaExGensXpath << mzr::eltName::model
		     << slash
		     << mzr::eltName::reactionGens
		     << slash
		     << eltName::gpaExchangeGen
		     << slash
		     << plx::eltName::plex;
      rPlexUnit.addOmniXpath(gpaExGensXpath.str());
    }

    void
    parseDomInput(xmlpp::Element* pRootElement,
		  xmlpp::Element* pModelElement,
		  xmlpp::Element* pStreamsElement,
		  xmlpp::Element* pEventsElement) throw(std::exception);

    void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
    {
      // This unit is just going to provide reaction generators,
      // and for now, I'm not emitting reaction generators in
      // a "state dump."
    }
  };
}

#endif
