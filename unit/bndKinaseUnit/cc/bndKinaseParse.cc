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

#include "mol/modMol.hh"
#include "bndKinase/modFam.hh"
#include "bndKinase/bndKinaseFam.hh"
#include "bndKinase/bndKinaseUnit.hh"
#include "bndKinase/bndKinaseEltName.hh"
#include "bndKinase/bndKinaseXcpt.hh"
#include "bndKinase/parseModGen.hh"
#include "bndKinase/parseBndKinaseGen.hh"
#include "bndKinase/parseBndOmniGen.hh"
#include "plex/plexDomParse.hh"

namespace bndKinase
{
  void
  bndKinaseUnit::parseDomInput(xmlpp::Element* pRootElement,
			       xmlpp::Element* pModelElement,
			       xmlpp::Element* pStreamsElement,
			       xmlpp::Element* pEventsElement)
    throw(std::exception)
  {
    // This unit only adds a reaction generator for now.
    xmlpp::Element* pReactionGensElt
      = domUtils::mustGetUniqueChild(pModelElement,
				     mzr::eltName::reactionGens);

    // Get the bndKinaseGen nodes.
    xmlpp::Node::NodeList bndKinaseGenNodes
      = pReactionGensElt->get_children(eltName::bndKinaseGen);

    // Parse bndKinase reaction generators.
    std::for_each(bndKinaseGenNodes.begin(),
		  bndKinaseGenNodes.end(),
		  parseBndKinaseGen(rMzrUnit,
				    rMolUnit,
				    rPlexUnit));

    // Get the modGen nodes.
    xmlpp::Node::NodeList modGenNodes
      = pReactionGensElt->get_children(eltName::modGen);

    // Add modGen reaction family for each of the generator nodes.
    std::for_each(modGenNodes.begin(),
		  modGenNodes.end(),
		  parseModGen(rMzrUnit,
			      rMolUnit,
			      rPlexUnit));

    // Get the bndOmniGen nodes.
    xmlpp::Node::NodeList bndOmniGenNodes
      = pReactionGensElt->get_children(eltName::bndOmniGen);

    // Add bndOmniFam reaction family for each of the generator nodes.
    std::for_each(bndOmniGenNodes.begin(),
		  bndOmniGenNodes.end(),
		  parseBndOmniGen(rMzrUnit,
				  rMolUnit,
				  rPlexUnit));
  }
}
