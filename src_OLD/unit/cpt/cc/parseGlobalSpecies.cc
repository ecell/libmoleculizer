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

#include "cpt/parseGlobalSpecies.hh"
#include "cpt/cptEltName.hh"

namespace cpt
{
  class parseBoundaryRate :
    public std::unary_function<xmlpp::Node*, void>
  {
    const compartmentGraph& rGraph;
    std::vector<double>& rBoundaryRates;

  public:
    parseBoundaryRate(const compartmentGraph& rCompartmentGraph,
		      std::vector<double>& rBoundaryRatesVector) :
      rGraph(rCompartmentGraph),
      rBoundaryRates(rBoundaryRatesVector)
    {}

    void
    operator()(xmlpp::Node* pBoundaryRateNode) const
      throw(utl::xcpt)
    {
      xmlpp::Element* pBoundaryRateElt
	= utl::dom::mustBeElementPtr(pBoundaryRateNode);

      std::string leftCompartmentName
	= utl::dom::mustGetAttrString(pBoundaryRateElt,
				      eltName::boundary_leftCompartmentAttr);
      compartment* pLeftCompartment
	= rGraph.mustFindCompartment(leftCompartmentName,
				     pBoundaryRateElt);
      
      std::string rightCompartmentName
	= utl::dom::mustGetAttrString(pBoundaryRateElt,
				      eltName::boundary_rightCompartmentAttr);
      compartment* pRightCompartment
	= rGraph.mustFindCompartment(rightCompartmentName,
				     pBoundaryRateElt);

      int boundaryIndex
	= rGraph.mustFindBoundaryIndex(pLeftCompartment,
				       pRightCompartment,
				       pBoundaryRateElt);

      double rate
	= utl::dom::mustGetAttrNNDouble(pBoundaryRateElt,
					eltName::boundaryRate_valueAttr);

      rBoundaryRates[boundaryIndex] = rate;
    }
  };

  std::vector<double>
  parseBoundaryRates(xmlpp::Node* pParentNode,
		     const compartmentGraph& rGraph)
  {
    xmlpp::Element* pParentElt
      = utl::dom::mustBeElementPtr(pParentNode);
      
    // Parse the default diffusion rate.
    xmlpp::Element* pDefaultDiffusionRateElt
      = utl::dom::mustGetUniqueChild(pParentElt,
				     eltName::defaultDiffusionRate);
    double defaultDiffusionRate
      = utl::dom::mustGetAttrDouble(pDefaultDiffusionRateElt,
				    eltName::defaultDiffusionRate_valueAttr);

    // Parse "exceptional" diffusion rates across specified boundaries.
    std::vector<double> boundaryRates(rGraph.boundaries.size(),
				      defaultDiffusionRate);

    xmlpp::Node::NodeList boundaryRateNodes
      = pParentElt->get_children(eltName::boundaryRate);

    std::for_each(boundaryRateNodes.begin(),
		  boundaryRateNodes.end(),
		  parseBoundaryRate(rGraph,
				    boundaryRates));

    return boundaryRates;
  }

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
    operator()(xmlpp::Node* pCompartmentPopNode) const
    {
      xmlpp::Element* pCompartmentPopElt
	= utl::dom::mustBeElementPtr(pCompartmentPopNode);

      std::string compartmentName
	= utl::dom::mustGetAttrString
	(pCompartmentPopElt,
	 eltName::compartmentPop_compartmentNameAttr);

      int compartmentNdx
	= rGraph.mustFindCompartmentIndex(compartmentName,
					  pCompartmentPopElt);

      int compartmentPop
	= utl::dom::mustGetAttrInt(pCompartmentPopElt,
				   eltName::compartmentPop_countAttr);

      rPops[compartmentNdx] = compartmentPop;
    }
  };

  std::vector<int>
  parseCompartmentPops(xmlpp::Node* pParentNode,
		       const compartmentGraph& rGraph)
  {
    xmlpp::Element* pParentElt
      = utl::dom::mustBeElementPtr(pParentNode);
    
    // Parse the default compartment population.
    xmlpp::Element* pDefaultCompartmentPopElt
      = utl::dom::mustGetUniqueChild(pParentElt,
				     eltName::defaultCompartmentPop);
    int defaultCompartmentPop
      = utl::dom::mustGetAttrInt(pDefaultCompartmentPopElt,
				 eltName::defaultCompartmentPop_countAttr);

    // Parse "exceptional" compartment populations.
    std::vector<int> compartmentPops(rGraph.compartments.size(),
				     defaultCompartmentPop);

    xmlpp::Node::NodeList compartmentPopNodes
      = pParentElt->get_children(eltName::compartmentPop);

    std::for_each(compartmentPopNodes.begin(),
		  compartmentPopNodes.end(),
		  parseCompartmentPop(rGraph,
				      compartmentPops));
    return compartmentPops;
  }
}
