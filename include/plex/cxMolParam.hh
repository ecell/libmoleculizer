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

#ifndef CXMOLPARAM_H
#define CXMOLPARAM_H

#include "mzr/feature.hh"
#include "mol/molState.hh"
#include "plex/plexFamily.hh"

namespace plx
{
  // This class adds functionality to a mol in the context of a complex
  // species.
  class cxMol :
    public mzr::featureContext<plexSpecies, plexMolSpec>
  {
  public:
    cxMol(const mzr::featureContext<plexSpecies, plexMolSpec>& rfc) :
      mzr::featureContext<plexSpecies, plexMolSpec>(rfc)
    {}

    // Get the index of the mol (plexMolSpec) in the plexFamily.
    plexMolSpec
    getMolSpec(void) const
    {
      return getSpec();
    }

    // Gets the population of the complex in which the mol occurs.
    // Used in most propensity calculations.
    int
    getPop(void) const
    {
      // Same as the template function.
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

    // Get the plex family in which the mol occurs.
    plexFamily&
    getPlexFamily(void) const
    {
      // Same as the template function.
      return getSpecies()->getFamily();
    }

    // Get the full plexParam of the complex in which the mol occurs.
    // This is typically used to update a substrate plex species
    // one of whose members is going away.
    const plexParam&
    getPlexParam(void) const
    {
      return getSpecies()->getParam();
    }

    // Extracts the vector of molParams from the plexParam retrieved
    // above.  This typically used as the beginning of the construction
    // of the plexParam of a product complex by means of allostery.
    const std::vector<bnd::molParam>&
    getMolParams(void) const
    {
      return getPlexParam().molParams;
    }

    // Gets the mol param of the "focus" mol by extracting it from the
    // above vector.  Use this when the reaction needs to know the state
    // of the mol.  The mol is capable of converting a corresponding
    // molParam into a vector of siteParams describing allosteric
    // binding when the molecule is in the specified state.
    bnd::molParam
    getMolParam(void) const
    {
      return getMolParams()[getMolSpec()];
    }
  };
}

#endif
