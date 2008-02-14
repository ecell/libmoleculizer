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

#ifndef CXSITEPARAM_H
#define CXSITEPARAM_H

#include "fnd/featureContext.hh"
#include "cpx/molState.hh"
#include "cpx/ftrSpec.hh"
#include "cpx/siteToShapeMap.hh"

namespace cpx
{
  // This class describes a binding site in the context  of a plexSpecies by
  // giving the index of the mol and the index of the binding site on the mol,
  // i.e. a siteSpec.
  template<class plexSpeciesT, class plexFamilyT>
  class cxSite :
    public fnd::featureContext<plexSpeciesT, siteSpec>
  {
  public:
    typedef plexSpeciesT plexSpeciesType;
    typedef plexFamilyT plexFamilyType;
    
    cxSite(typename cxSite::plexSpeciesType* pPlexSpecies,
	   const siteSpec& rSpec);

    // Get the site spec of the (free) site.
    siteSpec
    getSiteSpec(void) const;

    // Used in almost all propensity calculations.
    int
    getPop(void) const;

    // Used in almost all propensity calculations.
    double
    getPlexWeight(void) const;

    // Get the plex family on whose members the site appears.
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

    // Get the complex's siteParam for the site.  What is the shape
    // of the site in the current context?
    siteParam
    getSiteParam(void) const;
  };
}

#include "cpx/cxSiteImpl.hh"

#endif
