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
#include "clx/cptPlexSpecies.hh"
#include "clx/cptPlexFamily.hh"
#include "clx/clxEltName.hh"

namespace clx
{
  cptPlexSpecies::
  cptPlexSpecies(cptPlexFamily& rContainingFamily,
		 const cpx::siteToShapeMap& rSiteParams,
		 const std::vector<cpx::molParam>& rMolParams,
		 const cpt::compartmentGraph& rCompartmentGraph,
		 const std::vector<double>& rDiffusionRates) :
    cpt::globalSpecies(rCompartmentGraph,
		       rDiffusionRates),
    cpx::plexSpeciesMixin<cptPlexFamily>(rContainingFamily,
					 rSiteParams,
					 rMolParams)
  {}

  double
  cptPlexSpecies::
  getWeight(void) const
  {
    return cpx::plexSpeciesMixin<cptPlexFamily>::getWeight();
  }

  void
  cptPlexSpecies::
  notify(int notifyDepth)
  {
    rFamily.respond
      (fnd::newSpeciesStimulus<cptPlexSpecies>(this,
					       notifyDepth));
  }
  
  void
  cptPlexSpecies::
  respond(const fnd::newSpeciesStimulus<cpt::compartmentSpecies>& rStim)
  {
    ensureNotified(rStim.getNotificationDepth());
  }

  std::string
  cptPlexSpecies::
  getTag(void) const
  {
    return utl::stringify<const cptPlexSpecies*>(this);
  }

  std::string
  cptPlexSpecies::
  getName(void) const
  {
    return getInformativeName();
  }

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

xmlpp::Element*
  cptPlexSpecies::
  insertElt(xmlpp::Element* pExplicitSpeciesElt) const
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
	cml::cptMol* pMol = rFamily.getParadigm().mols[molNdx];
      
	pMol->insertInstanceState(pInstanceStatesElt,
				  molNdx,
				  molParams[molNdx]);
      }

    xmlpp::Element* pPopulationElt
      = pTaggedPlexSpeciesElt->add_child(eltName::population);

    pPopulationElt->set_attribute(eltName::population_countAttr,
				  utl::stringify<int>(getTotalPop()));

    // Adding redundant concentration element for use by ODE solver.  An
    // alternative would be to convert population to concentration (using
    // Java?)  during translation of state dump for ODE solver.
    //
    // This all needs to be broken out into separate compartments.
    double concentration
      = getTotalPop()/getCompartmentGraph().getTotalVolume();

    xmlpp::Element* pConcentrationElt
      = pTaggedPlexSpeciesElt->add_child(eltName::concentration);
    pConcentrationElt->set_attribute(eltName::concentration_valueAttr,
				     utl::stringify<double>(concentration));

    // Insert element for each compartment species, by compartment name
    // as well as index, giving the population and concentration.
    std::for_each(compartmentSpeciesVector.begin(),
		  compartmentSpeciesVector.end(),
		  insertCompartmentPop(rGraph,
				       pTaggedPlexSpeciesElt));

    // Add the updated flag for use by parametrizer.
    if(hasNotified())
      {
	pTaggedPlexSpeciesElt->add_child(eltName::updated);
      }

    return pTaggedPlexSpeciesElt;
  }
}


