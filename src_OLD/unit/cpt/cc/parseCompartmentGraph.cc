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

#include "cpt/parseCompartmentGraph.hh"
#include "cpt/cptEltName.hh"
#include "cpt/unitsMgr.hh"
#include "cpt/cptUnit.hh"

namespace cpt
{
  class parseCompartment :
    public std::unary_function<xmlpp::Node*, void>
  {
    compartmentGraph& rGraph;
  public:
    parseCompartment(compartmentGraph& rCompartmentGraph) :
      rGraph(rCompartmentGraph)
    {}

    void
    operator()(xmlpp::Node* pCompartmentNode) const
      throw(utl::xcpt)
    {
      xmlpp::Element* pCompartmentElt
	= utl::dom::mustBeElementPtr(pCompartmentNode);

      std::string compartmentName
	= utl::dom::mustGetAttrString(pCompartmentElt,
				      eltName::compartment_nameAttr);
      double compartmentVolume
	= utl::dom::mustGetAttrDouble(pCompartmentElt,
				      eltName::compartment_volumeAttr);

      // The compartment graph manages both its compartments and the
      // boundaries between them.
      rGraph.mustAddCompartment(new compartment(compartmentName,
						compartmentVolume),
				pCompartmentNode);
    }
  };

  class parseBoundary :
    public std::unary_function<xmlpp::Node*, void>
  {
    compartmentGraph& rGraph;
  public:
    parseBoundary(compartmentGraph& rCompartmentGraph) :
      rGraph(rCompartmentGraph)
    {}

    void
    operator()(xmlpp::Node* pBoundaryNode) const
    {
      xmlpp::Element* pBoundaryElt
	= utl::dom::mustBeElementPtr(pBoundaryNode);

      std::string leftCompartmentName
	= utl::dom::mustGetAttrString(pBoundaryElt,
				      eltName::boundary_leftCompartmentAttr);
      compartment* pLeftCompartment
	= rGraph.mustFindCompartment(leftCompartmentName,
				     pBoundaryElt);
      
      std::string rightCompartmentName
	= utl::dom::mustGetAttrString(pBoundaryElt,
				      eltName::boundary_rightCompartmentAttr);
      compartment* pRightCompartment
	= rGraph.mustFindCompartment(rightCompartmentName,
				     pBoundaryElt);

      double area
	= utl::dom::mustGetAttrPosDouble(pBoundaryElt,
				      eltName::boundary_areaAttr);

      double thickness
	= utl::dom::mustGetAttrPosDouble(pBoundaryElt,
					 eltName::boundary_thicknessAttr);

      // The compartment graph manages both its compartments and its
      // boundaries.
      rGraph.mustAddBoundary(new boundary(pRightCompartment,
					  pLeftCompartment,
					  area,
					  thickness));
    }
  };
      
  void
  parseCompartmentGraph::
  operator()(xmlpp::Node* pCompartmentGraphNode) const
    throw(utl::xcpt)
  {
    compartmentGraph& rGraph
      = rApp.getUnits()->getCptUnit().getCompartmentGraph();
    
    xmlpp::Element* pCompartmentGraphElt
      = utl::dom::mustBeElementPtr(pCompartmentGraphNode);

    xmlpp::Element* pCompartmentsElt
      = utl::dom::mustGetUniqueChild(pCompartmentGraphElt,
				     eltName::compartments);

    xmlpp::Node::NodeList compartmentNodes
      = pCompartmentsElt->get_children(eltName::compartment);

    std::for_each(compartmentNodes.begin(),
		  compartmentNodes.end(),
		  parseCompartment(rGraph));

    xmlpp::Element* pBoundariesElt
      = utl::dom::mustGetUniqueChild(pCompartmentGraphElt,
				     eltName::boundaries);

    xmlpp::Node::NodeList boundaryNodes
      = pBoundariesElt->get_children(eltName::boundary);

    std::for_each(boundaryNodes.begin(),
		  boundaryNodes.end(),
		  parseBoundary(rGraph));

    
  }
}
