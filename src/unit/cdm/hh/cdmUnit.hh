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

#ifndef CDM_CDMUNIT_H
#define CDM_CDMUNIT_H

#include "utl/dom.hh"
#include "cpt/cptUnit.hh"
#include "cml/cmlUnit.hh"
#include "clx/clxUnit.hh"
#include "cdm/cdmEltName.hh"
#include "cdm/dimerizeRxnGen.hh"
#include "cdm/decompRxnGen.hh"

namespace cdm
{

  /*! \ingroup dimerGroup
    \brief Provides dimerization and decomposition reactions.

    This unit provides the basic reactions for binding complexes
    together at free binding sites, and conversely, for breaking complexes
    apart at bindings.

    This unit also provides databases of rates constants connected
    with bindings and with %pairs of compatible binding sites.  These
    databases are available through static functions.  They are used by
    plexFamily routines in addition to being used in this unit, so
    that for the time being, this is an essentially mandatory unit. */
  class cdmUnit :
    public cpt::unit
  {
  public:

    cpt::cptUnit& rCptUnit;
    cml::cmlUnit& rCmlUnit;
    clx::clxUnit& rClxUnit;


    cdmUnit(cpt::cptApp& rCptApp,
	    cpt::cptUnit& refCptUnit,
	    cml::cmlUnit& refCmlUnit,
	    clx::clxUnit& refClxUnit) :
      cpt::unit("dimer",
		rCptApp),
      rCptUnit(refCptUnit),
      rCmlUnit(refCmlUnit),
      rClxUnit(refClxUnit)
    {
      // This unit isn't responsible for any model elements.

      // Register that this unit is responsible for dimerization/decomposition
      // generators.
      inputCap.addReactionGenName(eltName::dimerizationGen);

      // This unit isn't responsible for any explicit species, species streams,
      // or events.
    }

    void
    parseDomInput(xmlpp::Element* pRootElement,
		  xmlpp::Element* pModelElement,
		  xmlpp::Element* pStreamsElement,
		  xmlpp::Element* pEventsElement) throw(std::exception);

    void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception);
  };
}

#endif // CDM_CDMUNIT_H
