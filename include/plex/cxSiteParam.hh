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

#include "mzr/feature.hh"
#include "mol/siteShape.hh"
#include "mol/molState.hh"
#include "plex/plexFamily.hh"

namespace plx
{
  class cxSite :
    public mzr::featureContext<plexSpecies, plexSiteSpec>
  {
  public:
    // Does the wrapping.  No new data.
    cxSite(const mzr::featureContext<plexSpecies, plexSiteSpec>& rfc) :
      mzr::featureContext<plexSpecies, plexSiteSpec>(rfc)
    {}

    // Get the site spec of the (free) site.
    plexSiteSpec
    getSiteSpec(void) const
    {
      // Same as the template function.
      return getSpec();
    }

    // Used in almost all propensity calculations.
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

    // Get the plex family on whose members the site appears.
    plexFamily&
    getPlexFamily(void) const
    {
      // Same as the template function.
      return getSpecies()->getFamily();
    }

    // Get the full plexParam of the complex on which the site occurs.
    const plexParam&
    getPlexParam(void) const
    {
      // Perhaps species should return a reference, too.
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

    // Get the complex's siteParam for the site.  What is the shape
    // of the site in the current context?
    bnd::siteParam
    getSiteParam(void) const
    {
      const std::map<plexSiteSpec, bnd::siteParam>& rSiteParams
	= getPlexParam().siteParams;
      
      std::map<plexSiteSpec, bnd::siteParam>::const_iterator iSpecParam
	= rSiteParams.find(getSiteSpec());

      // Debugging code.
      if(iSpecParam == rSiteParams.end()) throw badSiteSpecXcpt();
      
      return iSpecParam->second;
    }
  };
}

#endif
