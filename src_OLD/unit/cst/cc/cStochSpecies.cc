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

#include "utl/string.hh"
#include "cpt/cptEltName.hh"
#include "cpt/globalSpecies.hh"
#include "cst/cstEltName.hh"
#include "cst/cStochSpecies.hh"

namespace cst
{
  // Inserts the compartment populations of a globalSpecies into
  // state output.
  class insertCompartmentPop :
    public std::unary_function<const cpt::compartmentSpecies*, void>
  {
    const cpt::compartmentGraph& rGraph;
    xmlpp::Element* pParent;

  public:
    insertCompartmentPop(const cpt::compartmentGraph& rCompartmentGraph,
			 xmlpp::Element* pTaggedStochSpeciesElt) :
      rGraph(rCompartmentGraph),
      pParent(pTaggedStochSpeciesElt)
    {}

    void
    operator()(const cpt::compartmentSpecies* pCompartmentSpecies) const
    {
      // Add element for this compartment.
      xmlpp::Element* pCompartmentPopElt
	= pParent->add_child(eltName::compartmentPop);

      // Add compartment name as attribute.
      pCompartmentPopElt->
	set_attribute(eltName::compartmentPop_compartmentNameAttr,
		      pCompartmentSpecies->getCompartment()->getName());

      // Add compartment population as attribute.
      pCompartmentPopElt->
	set_attribute(eltName::compartmentPop_populationAttr,
		      utl::stringify<int>(pCompartmentSpecies->getPop()));

      // Add compartment concentration as attribute.
      pCompartmentPopElt->
	set_attribute(eltName::compartmentPop_concentrationAttr,
		      utl::stringify<double>(pCompartmentSpecies->getConc()));
    }
  };

  // Elements describing a globalSpecies into state output.
  xmlpp::Element*
  cStochSpecies::
  insertElt(xmlpp::Element* pTaggedSpeciesElt) const
    throw(std::exception)
  {
    xmlpp::Element* pTaggedStochSpeciesElt
      = pTaggedSpeciesElt->add_child(eltName::taggedStochSpecies);

    pTaggedStochSpeciesElt->set_attribute(eltName::taggedStochSpecies_tagAttr,
					  getTag());

    pTaggedStochSpeciesElt->set_attribute(eltName::taggedStochSpecies_nameAttr,
					  getName());

    xmlpp::Element* pWeightElt
      = pTaggedStochSpeciesElt->add_child(eltName::weight);

    pWeightElt->set_attribute(eltName::weight_daltonsAttr,
			      utl::stringify<double>(getWeight()));

    // Insert element for each compartment species, by compartment name
    // as well as index, giving the population and concentration.
    std::for_each(compartmentSpeciesVector.begin(),
		  compartmentSpeciesVector.end(),
		  insertCompartmentPop(rGraph,
				       pTaggedStochSpeciesElt));

    return pTaggedStochSpeciesElt;
  }
}
