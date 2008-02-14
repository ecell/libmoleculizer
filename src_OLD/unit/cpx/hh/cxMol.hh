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

#include "fnd/featureContext.hh"
#include "cpx/ftrSpec.hh"
#include "cpx/prm.hh"

namespace cpx
{
  // This class adds functionality to a mol in the context of a complex
  // species.
  template<class plexSpeciesT, class plexFamilyT>
  class cxMol :
    public fnd::featureContext<plexSpeciesT, molSpec>
  {
  public:
    typedef plexSpeciesT plexSpeciesType;
    typedef plexFamilyT plexFamilyType;

    cxMol(plexSpeciesType* pPlexSpecies,
	  const molSpec& rSpec);

    // Get the index of the mol (molSpec) in the plexFamily.
    molSpec
    getMolSpec(void) const;

    // Gets the population of the complex in which the mol occurs.
    // Used in most propensity calculations.
    int
    getPop(void) const;

    // Used in almost all propensity calculations.
    double
    getPlexWeight(void) const;

    // Get the plex family in which the mol occurs.
    plexFamilyT&
    getPlexFamily(void) const;

    // Extracts the site shapes from the plexSpecies.
    const siteToShapeMap&
    getSiteToShapeMap(void) const;

    // Extracts the vector of molParams from the plexParam retrieved
    // above.  This typically used as the beginning of the construction
    // of the plexParam of a product complex by means of allostery.
    const std::vector<molParam>&
    getMolParams(void) const;

    // Gets the mol param of the "focus" mol by extracting it from the
    // above vector.  Use this when the reaction needs to know the state
    // of the mol.  The mol is capable of converting a corresponding
    // molParam into a vector of siteParams describing allosteric
    // binding when the molecule is in the specified state.
    molParam
    getMolParam(void) const;
  };
}

#include "cpx/cxMolImpl.hh"

#endif
