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

#ifndef PLEXSPEC_H
#define PLEXSPEC_H

/*! \file plexSpec.hh
  \ingroup plexSpeciesGroup
  \ingroup plexFeatureGroup

  \brief Feature specifiers, bindings, and other plex components. */

#include <iostream>
#include <utility>
#include "domUtils/domUtils.hh"

namespace plx
{
  class plex;

  /*! \ingroup plexSpeciesGroup
    \ingroup plexFeatureGroup
    \brief Specifier for a site on a complex.

    The pair (m, n) specifies the nth site on the mth mol in the complex. */
  class plexSiteSpec :
    public std::pair<int, int>
  {
  public:
    plexSiteSpec(void) :
      std::pair<int, int>(-1, -1)
    {}
      
    plexSiteSpec(int moleculeIndex,
		 int siteIndex) :
      std::pair<int, int>(moleculeIndex, siteIndex)
    {}

    int molNdx(void) const
    {
      return first;
    }

    void setMolNdx(int moleculeIndex)
    {
      first = moleculeIndex;
    }

    int siteNdx(void) const
    {
      return second;
    }

    void setSiteNdx(int siteIndex)
    {
      second = siteIndex;
    }

    xmlpp::Element*
    insertElt(xmlpp::Element* pBindingElt,
	      const plex& rPlex) const throw(std::exception);
  };

  /*! \ingroup plexSpeciesGroup
    \brief A binding in a complex.

    Consists of two site specs, specifying the two sites in the complex
    that are bound together. */
  class plexBinding : public std::pair<plexSiteSpec, plexSiteSpec>
  {
  public:
    plexBinding(void) :
      std::pair<plexSiteSpec, plexSiteSpec>(plexSiteSpec(), plexSiteSpec())
    {}
  
    plexBinding(const plexSiteSpec& rLeftSiteSpec,
		const plexSiteSpec& rRightSiteSpec) :
      std::pair<plexSiteSpec, plexSiteSpec>(rLeftSiteSpec, rRightSiteSpec)
    {}

    const plexSiteSpec&
    leftSite(void) const
    {
      return first;
    }

    plexSiteSpec&
    leftSite(void)
    {
      return first;
    }

    const plexSiteSpec&
    rightSite(void) const
    {
      return second;
    }

    plexSiteSpec&
    rightSite(void)
    {
      return second;
    }

    xmlpp::Element*
    insertElt(xmlpp::Element* pPlexElt,
	      const plex& rPlex) const throw(std::exception);
  };

  /*! \ingroup plexSpeciesGroup
  
  \brief Specifier of a binding in a complex.  Just the binding's index. */
  typedef int plexBindingSpec;

  /*! \ingroup plexSpeciesGroup

  \brief Specifier of a mol in a complex.  Just the mol's index. */
  typedef int plexMolSpec;
}

#endif // PLEXSPEC_H
