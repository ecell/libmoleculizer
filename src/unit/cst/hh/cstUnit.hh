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

#ifndef CST_CSTUNIT_H
#define CST_CSTUNIT_H

#include "utl/autoCatalog.hh"
#include "utl/dom.hh"
#include "cpt/cptUnit.hh"
#include "cpt/cptApp.hh"
#include "cst/cStochSpecies.hh"

namespace cst
{
  class cstUnit :
    public cpt::unit
  {
    utl::autoCatalog<cStochSpecies> stochSpeciesByName;

    // Can this go away?
    std::vector<cpt::cptReaction*> noSubstrateArrows;
  public:

    bool
    addStochSpecies(cStochSpecies* pStochSpecies);

    void
    mustAddStochSpecies(cStochSpecies* pStochSpecies,
			xmlpp::Node* pRequestingNode = 0)
      throw(utl::xcpt);
  
    cStochSpecies*
    findStochSpecies(const std::string& rSpeciesName) const;

    cStochSpecies*
    mustFindStochSpecies(const std::string& rSpeciesName,
			 xmlpp::Node* pRequestingNode = 0) const
      throw(utl::xcpt);

    void
    addNoSubstrateArrow(cpt::cptReaction* pArrow);

    cpt::cptUnit& rCptUnit;
    
    cstUnit(cpt::cptApp& rCptApp,
	    cpt::cptUnit& refCptUnit);
  
    void
    parseDomInput(xmlpp::Element* pRootElement,
		  xmlpp::Element* pModelElement,
		  xmlpp::Element* pStreamsElement,
		  xmlpp::Element* pEventsElement) throw(std::exception);

    void
    prepareToRun(xmlpp::Element* pRootElt,
		 xmlpp::Element* pModelElt,
		 xmlpp::Element* pStreamsElt,
		 xmlpp::Element* pEventsElt) throw(std::exception);

    void
    prepareToContinue(xmlpp::Element* pRootElt,
		      xmlpp::Element* pModelElt,
		      xmlpp::Element* pStreamsElt,
		      xmlpp::Element* pEventsElt,
		      std::map<std::string, std::string>& rTagToName,
		      xmlpp::Element* pTaggedSpeciesElement)
      throw(std::exception);

    void
    insertStateElts(xmlpp::Element* pRootElt)
      throw(std::exception);
  };
}

#endif // CST_CSTUNIT_H
