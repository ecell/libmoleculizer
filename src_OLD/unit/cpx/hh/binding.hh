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

#ifndef CPX_BINDING_H
#define CPX_BINDING_H

#include "cpx/ftrSpec.hh"

namespace cpx
{
  /*! \ingroup plexSpeciesGroup
    \brief A binding in a complex.

    Consists of two site specs, specifying the two sites in the complex
    that are bound together. */
  class binding : 
    public std::pair<siteSpec, siteSpec>
  {
  public:
    // I think that this (apparently useless, since there is no way to
    // modify the constructed binding) constructor was in order to
    // be able to make an std::vector<binding>.
    binding(void) :
      std::pair<siteSpec, siteSpec>(siteSpec(), siteSpec())
    {}
  
    binding(const siteSpec& rLeftSiteSpec,
	    const siteSpec& rRightSiteSpec) :
      std::pair<siteSpec, siteSpec>(rLeftSiteSpec, rRightSiteSpec)
    {}

    const siteSpec&
    leftSite(void) const
    {
      return first;
    }

    siteSpec&
    leftSite(void)
    {
      return first;
    }

    const siteSpec&
    rightSite(void) const
    {
      return second;
    }

    siteSpec&
    rightSite(void)
    {
      return second;
    }
  };
}

#endif // CPX_BINDING_H
