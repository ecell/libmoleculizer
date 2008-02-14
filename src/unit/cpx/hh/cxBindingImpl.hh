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

#ifndef CPX_CXBINDINGIMPL_H
#define CPX_CXBINDINGIMPL_H

namespace cpx
{
  template<class plexSpeciesT, class plexFamilyT>
  cxBinding<plexSpeciesT, plexFamilyT>::
  cxBinding(plexSpeciesType* pPlexSpecies,
	    const bindingSpec& rSpec) :
    fnd::featureContext<plexSpeciesType, bindingSpec>(pPlexSpecies,
						      rSpec)
  {}

  template<class plexSpeciesT, class plexFamilyT>
  bindingSpec
  cxBinding<plexSpeciesT, plexFamilyT>::
  getBindingSpec(void) const
  {
    return this->getSpec();
  }

  template<class plexSpeciesT, class plexFamilyT>
  int
  cxBinding<plexSpeciesT, plexFamilyT>::
  getPop(void) const
  {
    return this->getSpecies()->getPop();
  }

  template<class plexSpeciesT, class plexFamilyT>
  double
  cxBinding<plexSpeciesT, plexFamilyT>::
  getPlexWeight(void) const
  {
    return this->getSpecies()->getWeight();
  }

  template<class plexSpeciesT, class plexFamilyT>
  plexFamilyT&
  cxBinding<plexSpeciesT, plexFamilyT>::
  getPlexFamily(void) const
  {
    return this->getSpecies()->rFamily;
  }

  
  template<class plexSpeciesT, class plexFamilyT>
  const siteToShapeMap&
  cxBinding<plexSpeciesT, plexFamilyT>::
  getSiteToShapeMap(void) const
  {
    return this->getSpecies()->siteParams;
  }

  template<class plexSpeciesT, class plexFamilyT>
  const std::vector<molParam>&
  cxBinding<plexSpeciesT, plexFamilyT>::
  getMolParams(void) const
  {
    return this->getSpecies()->molParams;
  }

  template<class plexSpeciesT, class plexFamilyT>
  std::pair<siteParam, siteParam>
  cxBinding<plexSpeciesT, plexFamilyT>::
  getSiteParams(void) const
  {
    // Get the featured binding.
    binding theBinding
      = getPlexFamily().getParadigm().bindings[getBindingSpec()];

    // Get the shapes of the sites comprising the featured
    // binding in the featured plex species.
    siteParam leftSiteParam
      = getSiteToShapeMap().mustGetSiteShape(theBinding.leftSite());
    siteParam rightSiteParam
      = getSiteToShapeMap().mustGetSiteShape(theBinding.rightSite());

    return std::make_pair(leftSiteParam,
			  rightSiteParam);
  }
}

#endif // CPX_CXBINDINGIMPL_H
