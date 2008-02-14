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

#include "utl/dom.hh"
#include "utl/string.hh"
#include "plex/mzrPlexSpecies.hh"
#include "plex/mzrPlexFamily.hh"
#include "plex/plexEltName.hh"

namespace plx
{
  double
  mzrPlexSpecies::
  getWeight(void) const
  {
    return cpx::plexSpeciesMixin<mzrPlexFamily>::getWeight();
  }
  
  void
  mzrPlexSpecies::
  notify(int generateDepth)
  {
    rFamily.respond(fnd::newSpeciesStimulus<mzrPlexSpecies>(this,
							    generateDepth));
  }

  std::string
  mzrPlexSpecies::
  getName(void) const
  {
    return getInformativeName();
  }

  xmlpp::Element*
  mzrPlexSpecies::
  insertElt(xmlpp::Element* pExplicitSpeciesElt,
	    double molarFactor) const
    throw(std::exception)
  {
    // Insert tagged-plex-species element.
    xmlpp::Element* pTaggedPlexSpeciesElt
      = pExplicitSpeciesElt->add_child(eltName::taggedPlexSpecies);

    pTaggedPlexSpeciesElt->set_attribute(eltName::taggedPlexSpecies_tagAttr,
					 getTag());

    pTaggedPlexSpeciesElt->set_attribute(eltName::taggedPlexSpecies_nameAttr,
					 getName());

    // Insert the paradigm plex.
    rFamily.getParadigm().insertElt(pTaggedPlexSpeciesElt);

    // Insert the non-default instance states.
    //
    // Might be easier and non-harmful to insert all instance states.
    xmlpp::Element* pInstanceStatesElt
      = pTaggedPlexSpeciesElt->add_child(eltName::instanceStates);

    // I need molNdx to generate a pseudo instance name.
    for(int molNdx = 0;
	molNdx < (int) molParams.size();
	++molNdx)
      {
	bnd::mzrMol* pMol = rFamily.getParadigm().mols[molNdx];
      
	pMol->insertInstanceState(pInstanceStatesElt,
				  molNdx,
				  molParams[molNdx]);
      }

    xmlpp::Element* pPopulationElt
      = pTaggedPlexSpeciesElt->add_child(eltName::population);

    pPopulationElt->set_attribute(eltName::population_countAttr,
				  utl::stringify<int>(getPop()));

    // Adding redundant concentration element for use by ODE solver.  An
    // alternative would be to convert population to concentration (using
    // Java?)  during translation of state dump for ODE solver.
    double concentration = getPop()/molarFactor;

    xmlpp::Element* pConcentrationElt
      = pTaggedPlexSpeciesElt->add_child(eltName::concentration);
    pConcentrationElt->set_attribute(eltName::concentration_valueAttr,
				     utl::stringify<double>(concentration));

    // Add the updated flag for use by parametrizer.
    if(hasNotified())
      {
	pTaggedPlexSpeciesElt->add_child(eltName::updated);
      }

    return pTaggedPlexSpeciesElt;
  }
}
