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

#ifndef DIMERUNIT_H
#define DIMERUNIT_H

/*! \defgroup dimerGroup The dimer unit.
  \ingroup unitsGroup
  \brief Provides dimerization and decomposition reactions. */

/*! \file dimerUnit.hh
  \ingroup dimerGroup
  \brief Defines dimerUnit. */

#include "mzr/mzrUnit.hh"
#include "mol/molUnit.hh"
#include "plex/plexFeature.hh"
#include "dimer/dimerEltName.hh"
#include "dimer/dimerizeRxnGen.hh"
#include "dimer/decompRxnGen.hh"

namespace dimer
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
  class dimerUnit : public mzr::unit
  {
    // Reaction generators.
    mzr::autoVector<decompRxnGen> decompGens;
    mzr::autoVector<dimerizeRxnGenPair> dimerizeGens;

  public:

    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;


    void
    addDecompGen(decompRxnGen* pGenerator)
    {
      decompGens.addEntry(pGenerator);
    }

    void
    addDimerizeGen(dimerizeRxnGenPair* pGenerator)
    {
      dimerizeGens.addEntry(pGenerator);
    }

    dimerUnit(mzr::moleculizer& rMoleculizer,
	      mzr::mzrUnit& refMzrUnit,
	      bnd::molUnit& refMolUnit,
	      plx::plexUnit& refPlexUnit) :
      mzr::unit("dimer",
		rMoleculizer),
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
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

#endif
