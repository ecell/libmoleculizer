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
#include "cpt/globalReaction.hh"
#include "cpt/cptReaction.hh"

namespace cpt
{
  class addReactantToCpt :
    public std::unary_function<globalReaction::multMap::value_type, void>
  {
    cptReaction* pCptReaction;
  public:
    addReactantToCpt(cptReaction* pCompartmentReaction) :
      pCptReaction(pCompartmentReaction)
    {}

    void
    operator()(const argument_type& rGlobalSpeciesMultPair) const
    {
      int cptNdx = pCptReaction->getCompartmentIndex();

      globalSpecies* pGlobalReactant = rGlobalSpeciesMultPair.first;

      int multiplicity = rGlobalSpeciesMultPair.second;
      
      pCptReaction->addReactant(pGlobalReactant->getCompartmentSpecies(cptNdx),
				multiplicity);
    }
  };
    
  class addProductToCpt :
    public std::unary_function<globalReaction::multMap::value_type, void>
  {
    cptReaction* pCptReaction;
  public:
    addProductToCpt(cptReaction* pCompartmentReaction) :
      pCptReaction(pCompartmentReaction)
    {}

    void
    operator()(const argument_type& rGlobalSpeciesMultPair) const
    {
      int cptNdx = pCptReaction->getCompartmentIndex();

      globalSpecies* pGlobalProduct = rGlobalSpeciesMultPair.first;

      int multiplicity = rGlobalSpeciesMultPair.second;
      
      pCptReaction->addProduct(pGlobalProduct->getCompartmentSpecies(cptNdx),
				multiplicity);
    }
  };
    
  void
  globalReaction::
  finalizeCompartments(propensityDistro& rPropensityDistro)
  {
    int compartmentCount = rGraph.compartments.size();
    for(int compartmentNdx = 0;
	compartmentNdx < compartmentCount;
	++compartmentNdx)
      {
	// Construct a compartment reaction for this compartment.
	// These reactions, for the time being, are memory-managed
	// by the propensity distribution.
	cptReaction* pCptReaction
	  = new cptReaction(*this,
			    compartmentNdx);

	// Add reactants and products to the compartment reaction.
	//
	// Now, this routine adds this reaction to the sensitivity list
	// of each of its reactants as they are added, obsoleting
	// sensitizeToReactants
	std::for_each(reactants.begin(),
		      reactants.end(),
		      addReactantToCpt(pCptReaction));

	std::for_each(products.begin(),
		      products.end(),
		      addProductToCpt(pCptReaction));

	// Set the compartment reaction's rate.
	pCptReaction->setRate(rate);

	// Note that the function of "sensitizeToReactants" has now been
	// incorporated into "addReactant."

	// Add the reaction to the propensity distribution, but with
	// propensity 0.
	rPropensityDistro.addReaction(pCptReaction);

	// Update the reaction's propensity internally, in the
	// distribution, and update the totalPropensity.
	pCptReaction->respond(reactionStim(rPropensityDistro));
      }
  }

  void
  globalReaction::
  addReactant(globalSpecies* pSpecies,
	      int multiplicity)
  {
    fnd::basicReaction<globalSpecies>::addReactant(pSpecies,
						   multiplicity);
  }

  // State dump output for globalReactions.

  class insertReactantSpeciesRef :
    public std::unary_function<std::map<globalSpecies*, int>::value_type, void>
  {
    xmlpp::Element* pReactionElt;

  public:
    insertReactantSpeciesRef(xmlpp::Element* pReactionElement) :
      pReactionElt(pReactionElement)
    {}

    void
    operator()(const argument_type& rSpeciesMultPair) const
      throw(std::exception)
    {
      // Insert element for reaction reactant species.
      xmlpp::Element* pReactantSpeciesRefElt
	= pReactionElt->add_child(eltName::taggedSubstrate);

      // Add the name or tag of the substrate species as attribute.
      pReactantSpeciesRefElt
	->set_attribute(eltName::taggedSubstrate_tagAttr,
			rSpeciesMultPair.first->getTag());

      // Add multiplicity of substrate as attribute.
      pReactantSpeciesRefElt
	->set_attribute(eltName::taggedSubstrate_multAttr,
			utl::stringify<int>(rSpeciesMultPair.second));
    }
  };

  // It looks like I'm going back to emitting reaction substrates and
  // products, rather than substrates and deltas.  This is mainly for the
  // "immovable object," SBML.
  class insertProductSpeciesRef : public
  std::unary_function<std::map<globalSpecies*, int>::value_type, void>
  {
    xmlpp::Element* pReactionElt;

  public:
    insertProductSpeciesRef(xmlpp::Element* pReactionElement) :
      pReactionElt(pReactionElement)
    {}

    void
    operator()(const argument_type& rSpeciesMultPair) const
      throw(std::exception)
    {
      // Insert element for reaction product.
      xmlpp::Element* pProductSpeciesRefElt
	= pReactionElt->add_child(eltName::taggedProduct);

      // Add species tag as attribute.
      pProductSpeciesRefElt
	->set_attribute(eltName::taggedProduct_tagAttr,
			rSpeciesMultPair.first->getTag());

      // Add multiplicity as attribute.
      pProductSpeciesRefElt
	->set_attribute(eltName::taggedProduct_multAttr,
			utl::stringify<int>(rSpeciesMultPair.second));
    }
  };

  xmlpp::Element*
  globalReaction::
  insertElt(xmlpp::Element* pParentElt) const 
    throw(std::exception)
  {
    xmlpp::Element* pReactionElt
      = pParentElt->add_child(eltName::tagReaction);
      
    std::for_each(reactants.begin(),
		  reactants.end(),
		  insertReactantSpeciesRef(pReactionElt));

    std::for_each(products.begin(),
		  products.end(),
		  insertProductSpeciesRef(pReactionElt));

    // Additional scientific notation here for use by ECell.
    utl::dom::addDoubleParamChild(pReactionElt,
				  eltName::rate,
				  eltName::rate_valueAttr,
				  getRate());
    return pReactionElt;
  }

  int globalReaction::globalReactionCount = 0;
}
