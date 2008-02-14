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
#include "mzr/mzrEltName.hh"
#include "mol/molEltName.hh"
#include "stoch/stochEltName.hh"
#include "stoch/stochSpecies.hh"

namespace stoch
{
  xmlpp::Element*
  stochSpecies::insertElt(xmlpp::Element* pTaggedSpeciesElt,
			  double molarFactor) const
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

    xmlpp::Element* pPopulationElt
      = pTaggedStochSpeciesElt->add_child(eltName::population);

    pPopulationElt->set_attribute(eltName::population_countAttr,
				  utl::stringify<int>(getPop()));

    // Adding redundant concentration element for use by ODE solver.
    // An alternative would be to convert population to concentration
    // (using Java?) during translation of state dump fro ODE solver.
    double concentration
      = getPop() / molarFactor;

    xmlpp::Element* pConcentrationElt
      = pTaggedStochSpeciesElt->add_child(eltName::concentration);

    pConcentrationElt->set_attribute(eltName::concentration_valueAttr,
				     utl::stringify<double>(concentration));

    return pTaggedStochSpeciesElt;
  }
}
