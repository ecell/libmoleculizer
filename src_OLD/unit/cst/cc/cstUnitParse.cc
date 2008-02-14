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
#include "cpt/parseGlobalSpecies.hh"

#include "cst/cStochSpecies.hh"
#include "cst/cstUnit.hh"
#include "cst/cstEltName.hh"
#include "cst/badStochSpeciesTagXcpt.hh"

namespace cst
{
  class parseStochSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    cstUnit& rCstUnit;
    const cpt::compartmentGraph& rGraph;
    
  public:
    parseStochSpecies(cstUnit& refCstUnit,
		      const cpt::compartmentGraph& rCompartmentGraph) :
      rCstUnit(refCstUnit),
      rGraph(rCompartmentGraph)
    {}
    
    void
    operator()(xmlpp::Node* pStochSpeciesNode) const throw(std::exception)
    {
      xmlpp::Element* pStochSpeciesElt
	= utl::dom::mustBeElementPtr(pStochSpeciesNode);

      // Get the name of the stochastirator species.
      std::string speciesName
	= utl::dom::mustGetAttrString
	(pStochSpeciesElt,
	 eltName::stochSpecies_nameAttr);

      // Get the molecular weight of the stochastirator species.
      xmlpp::Element* pWeightElt
	= utl::dom::mustGetUniqueChild(pStochSpeciesElt,
				       eltName::weight);
      double molWeight
	= utl::dom::mustGetAttrPosDouble(pWeightElt,
					 eltName::weight_daltonsAttr);

      // Use core globalSpecies parsing routine to get diffusion rates over
      // the boundaries.
      std::vector<double>
	boundaryRates(parseBoundaryRates(pStochSpeciesElt,
					 rGraph));
      

      // Construct the stochastirator species.  These are memory-managed by
      // the cstUnit.
      cStochSpecies* pSpecies = new cStochSpecies(rGraph,
						  boundaryRates,
						  speciesName,
						  molWeight);

      // Register for memory management, construct dumpable, etc.
      rCstUnit.mustAddStochSpecies(pSpecies,
				   pStochSpeciesElt);
    }
  };

  void
  cstUnit::parseDomInput(xmlpp::Element* pRootElement,
			 xmlpp::Element* pModelElement,
			 xmlpp::Element* pStreamsElement,
			 xmlpp::Element* pEventsElement) throw(std::exception)
  {
    // Parse all the stoch species as named species.
    xmlpp::Element* pExplicitSpeciesElt
      = utl::dom::mustGetUniqueChild(pModelElement,
				     cpt::eltName::explicitSpecies);

    xmlpp::Node::NodeList stochSpeciesNodes
      = pExplicitSpeciesElt->get_children(eltName::stochSpecies);

    std::for_each(stochSpeciesNodes.begin(),
		  stochSpeciesNodes.end(),
		  parseStochSpecies(*this,
				    rCptUnit.getCompartmentGraph()));
  }
}
