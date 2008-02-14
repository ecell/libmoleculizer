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

#include "fnd/featureContext.hh"
#include "cpx/molState.hh"
#include "cpx/prm.hh"
#include "cpx/binding.hh"

namespace cpx
{
  // This class describes a binding in the context of a plexSpecies by
  // giving the index of the binding, and "dressing" it with accessors
  // to the underlying species, etc.
  template<class plexSpeciesT, class plexFamilyT>
  class cxBinding :
    public fnd::featureContext<plexSpeciesT,
			       bindingSpec>
  {
  public:
    typedef plexSpeciesT plexSpeciesType;
    typedef plexFamilyT plexFamilyType;

    cxBinding(plexSpeciesType* pPlexSpecies,
	      const bindingSpec& rSpec);

    // Get the binding spec (index) of the binding.  This is used to
    // extract the binding's parameter from the plex parameter, for
    // example.
    bindingSpec
    getBindingSpec(void) const;

    // Gets the population of the precise complex in which the
    // binding occurs.  Used in most propensity calculations.
    int
    getPop(void) const;

    // Used in almost all propensity calculations.
    double
    getPlexWeight(void) const;

    // Get the plex family in whose members this binding appears.
    plexFamilyT&
    getPlexFamily(void) const;

    // Extracts the site shapes from the plexSpecies.
    const siteToShapeMap&
    getSiteToShapeMap(void) const;

    // Extracts the vector of molParams from the plexParam.
    //
    // This is typically used in the construction of a product complex of a
    // reaction. Typically the mols of a product complex are in the same
    // states as they were in the reactant complexes, and this gets the states
    // of the mols in the reactant complexes.
    const std::vector<molParam>&
    getMolParams(void) const;

    // Gets the pair of binding site shapes connected with the "focus"
    // binding.
    // 
    // This is used to look up decomposition rates, for example, in the
    // decomposeExtrap reaction rate extrapolator.
    std::pair<siteParam, siteParam>
    getSiteParams(void) const;
  };
}

#include "cpx/cxBindingImpl.hh"

#endif
