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
#include "modKinase/modKinaseUnit.hh"
#include "modKinase/modKinaseUnitParse.hh"

namespace kinase
{
  void
  modKinaseUnit::parseDomInput(xmlpp::Element* pRootElement,
			       xmlpp::Element* pModelElement,
			       xmlpp::Element* pStreamsElement,
			       xmlpp::Element* pEventsElement)
    throw(std::exception)
  {
    // Get the header node for the reaction generators.
    xmlpp::Element* pReactionGensElt
      = domUtils::mustGetUniqueChild(pModelElement,
				     mzr::eltName::reactionGens);

    // Get the nucleotide binding generator elements, and add a reaction
    // family/generator for each.
    xmlpp::Node::NodeList nucBindRxnGenNodes
      = pReactionGensElt->get_children(eltName::nucleotideBindGen);

    std::for_each(nucBindRxnGenNodes.begin(),
		  nucBindRxnGenNodes.end(),
		  addNucBindRxnFam(rMzrUnit,
				   rMolUnit,
				   rStochUnit));

    // Get the kinase generator elements, and add a reaction family/generator
    // for each.
    xmlpp::Node::NodeList kinaseRxnGenNodes
      = pReactionGensElt->get_children(eltName::kinaseGen);

    std::for_each(kinaseRxnGenNodes.begin(),
		  kinaseRxnGenNodes.end(),
		  addKinaseRxnFam(rMzrUnit,
				  rMolUnit));

    // Get the phosphatase generator elements, and add a reaction
    // family/generator for each.
    xmlpp::Node::NodeList ptaseRxnGenNodes
      = pReactionGensElt->get_children(eltName::ptaseGen);

    std::for_each(ptaseRxnGenNodes.begin(),
		  ptaseRxnGenNodes.end(),
		  addPtaseRxnFam(rMzrUnit,
				 rMolUnit,
				 rStochUnit));
  }
}
