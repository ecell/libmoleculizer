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

#ifndef CXBINDINGPARAM_H
#define CXBINDINGPARAM_H

#include "mzr/featureContext.hh"
#include "mol/molState.hh"
#include "plex/plexSpecies.hh"
#include "plex/plexFamily.hh"

namespace plx
{
  class cxBinding :
    public mzr::featureContext<plexSpecies, plexBindingSpec>
  {
  public:
    cxBinding(const mzr::featureContext<plexSpecies, plexBindingSpec>& rfc) :
      mzr::featureContext<plexSpecies, plexBindingSpec>(rfc)
    {}

    // Get the binding spec (index) of the binding.  This is used to
    // extract the binding's parameter from the plex parameter, for
    // example.
    plexBindingSpec
    getBindingSpec(void) const
    {
      return getSpec();
    }

    // Gets the population of the precise complex in which the
    // binding occurs.  Used in most propensity calculations.
    int
    getPop(void) const
    {
      return getSpecies()->getPop();
    }

    // Used in almost all propensity calculations.
    double
    getPlexWeight(void) const
    {
      return getSpecies()->getWeight();
    }

    // Used in sensitizeToSubstrates.
    void
    addSensitiveReaction(mzr::reaction* pReaction) const
    {
      getSpecies()->addSensitiveReaction(pReaction);
    }

    // Get the plex family in whose members this binding appears.
    plexFamily&
    getPlexFamily(void) const
    {
      return getSpecies()->getFamily();
    }

    // Get the full plexParam of the specific complex in which the
    // binding occurs.  This is used typically in decrementing a
    // substrate species when a reaction happens.
    const plexParam&
    getPlexParam(void) const
    {
      return getSpecies()->getParam();
    }

    // Extracts the vector of molParams from the plexParam.  This
    // is typically used as the beginning of the construction of
    // the plexParam of a product complex by means of allostery.
    const std::vector<bnd::molParam>&
    getMolParams(void) const
    {
      return getPlexParam().molParams;
    }

    // Gets the pair of binding site shapes connected with this binding.
    // This is used to look up decomposition rates, for example, in the
    // decomposeExtrap reaction rate extrapolator.
    std::pair<bnd::siteParam, bnd::siteParam>
    getSiteParams(void) const
    {
      // Get the featured binding.
      plexBinding binding
	= getPlexFamily().getParadigm().bindings[getBindingSpec()];

      // Get the shapes of the sites comprising the featured
      // binding in the featured plex species.
      bnd::siteParam leftSiteParam
	= getPlexParam().getSiteParam(binding.leftSite());
      
      bnd::siteParam rightSiteParam
	= getPlexParam().getSiteParam(binding.rightSite());

      return std::make_pair(leftSiteParam,
			    rightSiteParam);
    }

    // Get the complex's bindingParam for this binding.  We need
    // this to compute the propensity of the binding to decompose.
//     bindingParam
//     getBindingParam(void) const
//     {
//       return getPlexParam().bindingParams[getBindingSpec()];
//     }

  };
}

#endif
