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
#include "cpt/unitsMgr.hh"
#include "cpt/cptUnit.hh"
#include "cpt/createEvent.hh"
#include "cpt/parseCreateEvent.hh"
#include "cpt/parseGlobalSpecies.hh"

namespace cpt
{
  // This class is for parsing the populations of a given species
  // in the different compartments.
  class parseCompartmentPop :
    public std::unary_function<xmlpp::Node*, void>
  {
    const compartmentGraph& rGraph;
    std::vector<int>& rPops;

  public:
    parseCompartmentPop(const compartmentGraph& rCompartmentGraph,
			std::vector<int>& rCompartmentPops) :
      rGraph(rCompartmentGraph),
      rPops(rCompartmentPops)
    {}

    void
    operator()(xmlpp::Node* pCompartmentPopNode)
    {
      xmlpp::Element* pCompartmentPopElt
	= utl::dom::mustBeElementPtr(pCompartmentPopNode);

      std::string compartmentName
	= utl::dom::mustGetAttrString(pCompartmentPopElt,
				      eltName::compartmentPop_compartmentNameAttr);
      int compartmentIndex
	= rGraph.mustFindCompartmentIndex(compartmentName,
					  pCompartmentPopElt);
      
      int compartmentPop
	= utl::dom::mustGetAttrNNInt(pCompartmentPopElt,
				     eltName::compartmentPop_countAttr);

      rPops[compartmentIndex] = compartmentPop;
    }
  };
    
  void
  parseCreateEvent::
  operator()(xmlpp::Node* pCreateEventNode)
    throw(utl::xcpt)
  {
    cptUnit& rCptUnit = rApp.getUnits()->getCptUnit();
    const compartmentGraph& rGraph = rCptUnit.getCompartmentGraph();
    
    xmlpp::Element* pCreateEventElt
      = utl::dom::mustBeElementPtr(pCreateEventNode);

    xmlpp::Element* pSpeciesRefElt
      = utl::dom::mustGetUniqueChild(pCreateEventElt,
				     eltName::speciesRef);
    std::string speciesName
      = utl::dom::mustGetAttrString(pSpeciesRefElt,
				    eltName::speciesRef_nameAttr);

    globalSpecies* pSpecies
      = rCptUnit.mustFindSpecies(speciesName,
				 pSpeciesRefElt);

    // Get the numbers of molecules to create in the different compartments.
    std::vector<int> compartmentPops
      = parseCompartmentPops(pCreateEventNode,
			     rGraph);

    // Construct the event and install it in its deletion pit.
    createEvent* pCreateEvent = new createEvent(pSpecies,
						compartmentPops);
    rCptUnit.addUserEvent(pCreateEvent);

    // Get the time at which to schedule the event.
    // 
    // Note that we will probably want to be able to schedule events
    // at negative times.
    double time
      = utl::dom::mustGetAttrDouble(pCreateEventElt,
				    eltName::createEvent_timeAttr);

    // Schedule the event if it is supposed to occur in the future.
    // (which might not be the case if this simulation is being
    // restarted from a state dump.)
    double now = rApp.getSimTime();
    if(now <= time)
      {
	rApp.getQueue().scheduleEvent(pCreateEvent,
				    time);
      }
  }
}
