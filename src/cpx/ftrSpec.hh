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

#ifndef CPX_FTRSPEC_H
#define CPX_FTRSPEC_H

#include <utility>

namespace cpx
{
  /*! \ingroup plexSpeciesGroup
    \ingroup plexFeatureGroup
    \brief Specifier for a site on a complex.

    The pair (m, n) specifies the nth site on the mth mol in the complex. */
  class siteSpec :
    public std::pair<int, int>
  {
  public:
    siteSpec(void) :
      std::pair<int, int>(-1, 
			  -1)
    {}
      
    siteSpec(int moleculeIndex,
		 int bindingSiteIndex) :
      std::pair<int, int>(moleculeIndex, 
			  bindingSiteIndex)
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
  };

  /*! \ingroup plexSpeciesGroup
  
  \brief Specifier of a binding in a complex.  Just the binding's index. */
  typedef int bindingSpec;

  /*! \ingroup plexSpeciesGroup

  \brief Specifier of a mol in a complex.  Just the mol's index. */
  typedef int molSpec;

}

#endif // CPX_FTRSPEC_H
