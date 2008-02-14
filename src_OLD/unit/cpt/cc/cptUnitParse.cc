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

#include "cpt/cptEltName.hh"
#include "cpt/parseCompartmentGraph.hh"
#include "cpt/parseGlobalReaction.hh"
#include "cpt/parseDumpStream.hh"
#include "cpt/parseCreateEvent.hh"
#include "cpt/parseStopEvent.hh"
#include "cpt/parseDumpStateEvent.hh"
#include "cpt/noStopEventWarning.hh"
#include "cpt/cptUnit.hh"

namespace cpt
{
  void
  cptUnit::
  parseDomInput(xmlpp::Element* pRootElement,
		xmlpp::Element* pModelElement,
		xmlpp::Element* pStreamsElement,
		xmlpp::Element* pEventsElement)
    throw(std::exception)
  {
    // Lots of stuff in all units appears to depend on the compartment graph
    // so (at least temporarily) moving parsing of the compartment graph
    // to the application's constructor
//     xmlpp::Element* pCompartmentGraphElt
//       = utl::dom::mustGetUniqueChild(pModelElement,
// 				     eltName::compartmentGraph);

//     // Parse compartment graph.
//     parseCompartmentGraph graphParser(rCptApp);
//     graphParser(pCompartmentGraphElt);

    // Parse explicit reactions.
    //
    // Eventually, there could be more than one kind of reaction,
    // so that this header element for all the different kinds isn't
    // useless.
    xmlpp::Element* pExplicitReactionsElt
      = utl::dom::mustGetUniqueChild(pModelElement,
				     eltName::explicitReactions);
    xmlpp::Node::NodeList reactionNodes
      = pExplicitReactionsElt->get_children(eltName::reaction);

    std::for_each(reactionNodes.begin(),
		  reactionNodes.end(),
		  parseGlobalReaction(rCptApp));

    // Parse dump-streams.
    //
    // At this point, there is only one kind of dump-stream, so the
    // dump-streams element looks redundant.
    xmlpp::Element* pDumpStreamsElt
      = utl::dom::mustGetUniqueChild(pStreamsElement,
				     eltName::dumpStreams);

    xmlpp::Node::NodeList dumpStreamNodes
      = pDumpStreamsElt->get_children(eltName::dumpStream);

    // Parse the dumpStream nodes to get the vector of tabDumpEvents, which
    // are special among events in being needed when state is dumped.
    std::for_each(dumpStreamNodes.begin(),
		  dumpStreamNodes.end(),
		  parseDumpStream(*this));

    // Parse create events.
    xmlpp::Node::NodeList createEventNodes
      = pEventsElement->get_children(eltName::createEvent);
    std::for_each(createEventNodes.begin(),
		  createEventNodes.end(),
		  parseCreateEvent(rCptApp));

    // Get stop event nodes.
    xmlpp::Node::NodeList stopEventNodes
      = pEventsElement->get_children(eltName::stopEvent);

    // Issue a warning if there are no stop event nodes.
    if(0 == stopEventNodes.size())
      noStopEventWarning(pEventsElement).warn();

    // Parse the stop events.
    std::for_each(stopEventNodes.begin(),
		  stopEventNodes.end(),
		  parseStopEvent(rCptApp));

    // Parse the dump-state events.
    // Parse dump-state-events.
    xmlpp::Node::NodeList dumpStateEventNodes
      = pEventsElement->get_children(eltName::dumpStateEvent);
    std::for_each(dumpStateEventNodes.begin(),
		  dumpStateEventNodes.end(),
		  parseDumpStateEvent(rCptApp));
  }
}
