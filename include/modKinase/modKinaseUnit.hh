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

#ifndef MODKINASEUNIT_H
#define MODKINASEUNIT_H

#include "mzr/mzrUnit.hh"
#include "mol/molUnit.hh"
#include "plex/plexUnit.hh"
#include "stoch/stochUnit.hh"
#include "modKinase/kinaseEltName.hh"

namespace kinase
{
  class modKinaseUnit : public mzr::unit
  {
  public:
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    stoch::stochUnit& rStochUnit;
    
    modKinaseUnit(mzr::moleculizer& rMoleculizer,
		  mzr::mzrUnit& refMzrUnit,
		  bnd::molUnit& refMolUnit,
		  plx::plexUnit& refPlexUnit,
		  stoch::stochUnit& refStochUnit) :
      mzr::unit("modKinase",
		rMoleculizer),
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      rStochUnit(refStochUnit)
    {
      inputCap.addReactionGenName(eltName::nucleotideBindGen);
      inputCap.addReactionGenName(eltName::kinaseGen);
      inputCap.addReactionGenName(eltName::ptaseGen);
    }

    // There are no omniplexes that need processing by the plex
    // unit.

    void
    parseDomInput(xmlpp::Element* pRootElement,
		  xmlpp::Element* pModelElement,
		  xmlpp::Element* pStreamsElement,
		  xmlpp::Element* pEventsElement) throw(std::exception);

    void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
    {
      // Since this unit only provides reaction generators, it has
      // no state to dump.
    }
  };
}

#endif // MODKINASEUNIT_H
