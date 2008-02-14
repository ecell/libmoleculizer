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
#include "cpt/parseGlobalReaction.hh"

namespace cpt
{
  class addProductGlobalSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    cptUnit& rCptUnit;
    globalReaction* pReaction;

  public:
    addProductGlobalSpecies(cptUnit& refCptUnit,
			    globalReaction* pParsedReaction) :
      rCptUnit(refCptUnit),
      pReaction(pParsedReaction)
    {}

    void
    operator()(xmlpp::Node* pProductGlobalSpeciesRefNode) const
      throw(std::exception)
    {
      xmlpp::Element* pProductGlobalSpeciesRefElt
	= utl::dom::mustBeElementPtr(pProductGlobalSpeciesRefNode);

      // Get the name of the substrate species.
      std::string globalSpeciesName
	= utl::dom::mustGetAttrString(pProductGlobalSpeciesRefElt,
				      eltName::productGlobalSpeciesRef_nameAttr);

      globalSpecies* pSpecies
	= rCptUnit.mustFindSpecies(globalSpeciesName,
				   pProductGlobalSpeciesRefElt);

      // Get the multiplicity of the product, which can be positive or
      // negative.
      int multiplicity
	= utl::dom::mustGetAttrPosInt(pProductGlobalSpeciesRefElt,
				      eltName::productGlobalSpeciesRef_multAttr);

      // Add the product species/multiplicity to the reaction.
      pReaction->addProduct(pSpecies,
			    multiplicity);
    }
  };

  class addReactantGlobalSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    cptUnit& rCptUnit;
    globalReaction* pReaction;

  public:
    addReactantGlobalSpecies(cptUnit& refCptUnit,
			     globalReaction* pParsedReaction) :
      rCptUnit(refCptUnit),
      pReaction(pParsedReaction)
    {}

    void
    operator()(xmlpp::Node* pReactantGlobalSpeciesRefNode) const
      throw(std::exception)
    {
      xmlpp::Element* pReactantGlobalSpeciesRefElt
	= utl::dom::mustBeElementPtr(pReactantGlobalSpeciesRefNode);

      // Get the name of the global reactant species.
      std::string globalSpeciesName
	= utl::dom::mustGetAttrString(pReactantGlobalSpeciesRefElt,
				      eltName::reactantGlobalSpeciesRef_nameAttr);

      globalSpecies* pSpecies
	= rCptUnit.mustFindSpecies(globalSpeciesName,
				   pReactantGlobalSpeciesRefElt);

      // Get the multiplicity of the reactant, which must be positive.
      int multiplicity
	= utl::dom::mustGetAttrPosInt(pReactantGlobalSpeciesRefElt,
				      eltName::reactantGlobalSpeciesRef_multAttr);

      // Add the reactant/multiplicity to the reaction.
      pReaction->addReactant(pSpecies,
			     multiplicity);
    }
  };

  void
  parseGlobalReaction::
  operator()(xmlpp::Node* pGlobalReactionNode) const
    throw(utl::xcpt)
  {
    cptUnit& rCptUnit = rApp.getUnits()->getCptUnit();

    // Construct the reaction, and register it for memory management.
    globalReaction* pParsedReaction
      = new globalReaction(rCptUnit.getCompartmentGraph());
    rCptUnit.addUserReaction(pParsedReaction);
    
    // Get the list of reactant nodes.
    xmlpp::Node::NodeList reactantGlobalSpeciesRefNodes
      = pGlobalReactionNode->get_children(eltName::reactantGlobalSpeciesRef);

    // Add the reactants to the reaction.
    std::for_each(reactantGlobalSpeciesRefNodes.begin(),
		  reactantGlobalSpeciesRefNodes.end(),
		  addReactantGlobalSpecies(rCptUnit,
					   pParsedReaction));

    // Get the list of product nodes.
    xmlpp::Node::NodeList productGlobalSpeciesRefNodes
      = pGlobalReactionNode->get_children(eltName::productGlobalSpeciesRef);

    // Add the products to the reaction.
    std::for_each(productGlobalSpeciesRefNodes.begin(),
		  productGlobalSpeciesRefNodes.end(),
		  addProductGlobalSpecies(rCptUnit,
					  pParsedReaction));

    // Get the reaction rate element.
    xmlpp::Element* pRateElt
      = utl::dom::mustGetUniqueChild(pGlobalReactionNode,
				     eltName::rate);
    double rate
      = utl::dom::mustGetAttrDouble(pRateElt,
				    eltName::rate_valueAttr);
    pParsedReaction->setRate(rate);

    // Make the compartment reactions and add them to the application's
    // propensity distribution.
    pParsedReaction->finalizeCompartments(rApp.getPropensities());
  }
}
