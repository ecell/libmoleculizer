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

#ifndef BNDKINASEUNIT_H
#define BNDKINASEUNIT_H

#include "mzr/mzrUnit.hh"
#include "mol/molUnit.hh"
#include "plex/plexUnit.hh"
#include "bndKinase/bndKinaseEltName.hh"

namespace bndKinase
{
  class bndKinaseUnit : public mzr::unit
  {
  public:
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;

    bndKinaseUnit(mzr::moleculizer& rMoleculizer,
		  mzr::mzrUnit& refMzrUnit,
		  bnd::molUnit& refMolUnit,
		  plx::plexUnit& refPlexUnit) :
      mzr::unit("bndKinase",
		rMoleculizer),
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {
      // Register reaction generator names.  This is used by the parser
      // to verify that all elements appearing in the input file are "claimed"
      // by some module.  Why did I think that was a good idea?  Now, "open
      // schema" i.e. pay no attention to extraneous matter while parsing.
      // seems more rational to me.
      inputCap.addReactionGenName(eltName::bndKinaseGen);
      inputCap.addReactionGenName(eltName::modGen);

      // Register the enabling complexes for bndKinase reaction gens as
      // omniplexes, for processing by the plexUnit.
      const std::string slash("/");
      std::ostringstream bndKinaseGensXpath;
      bndKinaseGensXpath << mzr::eltName::model
			 << slash
			 << mzr::eltName::reactionGens
			 << slash
			 << eltName::bndKinaseGen;
      rPlexUnit.addOmniXpath(bndKinaseGensXpath.str());

      // Register the enabling complexs for modification reaction gens as
      // omniplexes, for processing by the plexUnit.
      std::ostringstream modGensXpath;
      modGensXpath << mzr::eltName::model
		   << slash
		   << mzr::eltName::reactionGens
		   << slash
		   << eltName::modGen;
      rPlexUnit.addOmniXpath(modGensXpath.str());
    }

    // The input parsing routine for this unit.
    void
    parseDomInput(xmlpp::Element* pRootElement,
		  xmlpp::Element* pModelElement,
		  xmlpp::Element* pStreamsElement,
		  xmlpp::Element* pEventsElement)
      throw(std::exception);

    // The state output routine for this unit.
    void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
    {
      // For now, this unit just provides reaction generators, which
      // for now I'm not dumping in "state dump."
    }
  };
}

#endif // BNDKINASEUNIT_H
