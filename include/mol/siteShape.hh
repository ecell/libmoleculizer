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

#ifndef SITESHAPE_H
#define SITESHAPE_H

#include "domUtils/domUtils.hh"

namespace bnd
{
  /*! \ingroup plexSpeciesGroup

  \brief Placeholder for binding site shape.

  At this point in time, and possibly forever, dealing with the
  shape of a binding site in a realistic way might be impossible,
  but would certainly be slow.  So each site shape is just a unique
  thing with a name, and all binding propensities must be explicitly
  given by the user.

  Someday, binding site shapes could be shape-theoretic entities
  from which binding properties are actually deduced.  Someday,
  real protein geometry will have to be dealt with. */
  class siteShape
  {
  public:
    std::string name;

    siteShape(const std::string& rName) :
      name(rName)
    {}

    siteShape(xmlpp::Node* pSiteShapeNode) throw(std::exception);

    class constructor :
      public std::unary_function<xmlpp::Node*, siteShape>
    {
    public:
      siteShape
      operator()(xmlpp::Node* pSiteShapeNode) throw(std::exception)
      {
	return siteShape(pSiteShapeNode);
      }
    };

    const std::string&
    getName(void) const
    {
      return name;
    }
  };

  // This is a nouvelle regime change: siteShapes are generally const.
  typedef const siteShape* siteParam;
}

#endif // SITESHAPE_H
