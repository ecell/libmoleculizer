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

#ifndef CPX_OMNISTRUCTUREQUERYIMPL_H
#define CPX_OMNISTRUCTUREQUERYIMPL_H

#include "cpx/binding.hh"

namespace cpx
{
  class bindingContainsSite :
    public std::unary_function<binding, bool>
  {
    // The site being searched for among the bindings
    // AFTER it has been translated by the omni imbedding.
    const siteSpec& rQuerySiteSpec;

  public:
    bindingContainsSite(const siteSpec& rSiteSpec) :
      rQuerySiteSpec(rSiteSpec)
    {}
    
    bool
    operator()(const binding& rBinding) const
    {
      return (rBinding.leftSite() == rQuerySiteSpec)
	|| (rBinding.rightSite() == rQuerySiteSpec);
    }
  };

  template<class plexT>
  bool
  omniFreeSiteQuery<plexT>::
  operator()(const omniStructureQueryArg<plexT>& rArg) const
  {
    const plexT& rPlex = rArg.rPlex;
    const plexIso& rInjection = rArg.rInjection;

    // Translate the spec of the site we want to check into the indexing of
    // the given plex using the injection of the omniplex into the given plex
    // where it has been recognized.
    siteSpec translatedSiteSpec
      = rInjection.forward.applyToSiteSpec(freeSiteSpec);

    // Return true if no binding in the plex contains the translated site
    // spec; i.e. if the site specified by freeSiteSpec in the omni
    // is also free in the plex where the omni is embedded.
    return
      (rPlex.bindings.end()
       == std::find_if(rPlex.bindings.begin(),
		       rPlex.bindings.end(),
		       bindingContainsSite(translatedSiteSpec)));

  }
}

#endif // CPX_OMNISTRUCTUREQUERYIMPL_H
