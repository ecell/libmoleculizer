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

#ifndef SCAFFOLDUNIT_H
#define SCAFFOLDUNIT_H

#include "mzr/unit.hh"
#include "mzr/mzrEltName.hh"
#include "plex/plexUnit.hh"
#include "plex/plexEltName.hh"
#include "scaffold/scafEltName.hh"

namespace scaf
{
  class scaffoldUnit : public mzr::unit
  {
  public:
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
      
    scaffoldUnit(mzr::moleculizer& rMoleculizer,
		 mzr::mzrUnit& refMzrUnit,
		 bnd::molUnit& refMolUnit,
		 plx::plexUnit& refPlexUnit) :
      mzr::unit("scaffold",
		rMoleculizer),
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {
      // This unit provides only reacton generators, but all involve
      // omniplexes.
      inputCap.addReactionGenName(eltName::twentyElevenGen);
      inputCap.addReactionGenName(eltName::elevenSevenGen);
      inputCap.addReactionGenName(eltName::sevenThreeGen);
      inputCap.addReactionGenName(eltName::threeSevenGen);

      // Register omniplex nodes for processing by the plexUnit.
      const std::string slash("/");
      std::ostringstream twentyElevenGensXpath;
      twentyElevenGensXpath << mzr::eltName::model
			    << slash
			    << mzr::eltName::reactionGens
			    << slash
			    << eltName::twentyElevenGen
			    << slash
			    << plx::eltName::plex;
      rPlexUnit.addOmniXpath(twentyElevenGensXpath.str());

      std::ostringstream elevenSevenGensXpath;
      elevenSevenGensXpath << mzr::eltName::model
			   << slash
			   << mzr::eltName::reactionGens
			   << slash
			   << eltName::elevenSevenGen
			   << slash
			   << plx::eltName::plex;
      rPlexUnit.addOmniXpath(elevenSevenGensXpath.str());

      std::ostringstream sevenThreeGensXpath;
      sevenThreeGensXpath << mzr::eltName::model
			  << slash
			  << mzr::eltName::reactionGens
			  << slash
			  << eltName::sevenThreeGen
			  << slash
			  << plx::eltName::plex;
      rPlexUnit.addOmniXpath(sevenThreeGensXpath.str());

      std::ostringstream threeSevenGensXpath;
      threeSevenGensXpath << mzr::eltName::model
			  << slash
			  << mzr::eltName::reactionGens
			  << slash
			  << eltName::threeSevenGen
			  << slash
			  << plx::eltName::plex;
      rPlexUnit.addOmniXpath(threeSevenGensXpath.str());
    }

    void
    parseDomInput(xmlpp::Element* pRootElement,
		  xmlpp::Element* pModelElement,
		  xmlpp::Element* pStreamsElement,
		  xmlpp::Element* pEventsElement) throw(std::exception);

    void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
    {
      // For now, this unit just provides reaction generators, which
      // for now I'm not dumping in "state dump."
    }
  };
}

#endif // SCAFFOLDUNIT_H
