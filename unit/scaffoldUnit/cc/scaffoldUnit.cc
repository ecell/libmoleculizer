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

#include "mzr/mzrEltName.hh"
#include "scaffold/scaffoldUnit.hh"
#include "scaffold/scaffoldUnitParse.hh"
#include "scaffold/scafEltName.hh"

namespace scaf
{
  void
  scaffoldUnit::parseDomInput(xmlpp::Element* pRootElement,
			      xmlpp::Element* pModelElement,
			      xmlpp::Element* pStreamsElement,
			      xmlpp::Element* pEventsElement)
    throw(std::exception)
  {
    // Get the header node for all reaction generators.
    xmlpp::Element* pReactionGensElt
      = domUtils::mustGetUniqueChild(pModelElement,
				     mzr::eltName::reactionGens);

    // Get the twenty-eleven-gen nodes.
    xmlpp::Node::NodeList twentyElevenGenNodes
      = pReactionGensElt->get_children(eltName::twentyElevenGen);

    // Add twenty-eleven kinase reaction family for each.
    std::for_each(twentyElevenGenNodes.begin(),
		  twentyElevenGenNodes.end(),
		  addTwentyElevenRxnFam(rMzrUnit,
					rMolUnit,
					rPlexUnit));

    // Get the eleven-seven-gen nodes.
    xmlpp::Node::NodeList elevenSevenGenNodes
      = pReactionGensElt->get_children(eltName::elevenSevenGen);

    // Add eleven-seven kinase reaction family for each.
    std::for_each(elevenSevenGenNodes.begin(),
		  elevenSevenGenNodes.end(),
		  addElevenSevenRxnFam(rMzrUnit,
				       rMolUnit,
				       rPlexUnit));

    // Get the seven-three-gen nodes.
    xmlpp::Node::NodeList sevenThreeGenNodes
      = pReactionGensElt->get_children(eltName::sevenThreeGen);

    // Add seven-three kinase reaction family for each.
    std::for_each(sevenThreeGenNodes.begin(),
		  sevenThreeGenNodes.end(),
		  addSevenThreeRxnFam(rMzrUnit,
				      rMolUnit,
				      rPlexUnit));

    // Get the three-seven-gen nodes.
    xmlpp::Node::NodeList threeSevenGenNodes
      = pReactionGensElt->get_children(eltName::threeSevenGen);

    // Add three-seven kinase reaction family for each.
    std::for_each(threeSevenGenNodes.begin(),
		  threeSevenGenNodes.end(),
		  addThreeSevenRxnFam(rMzrUnit,
				      rMolUnit,
				      rPlexUnit));
  }
}
