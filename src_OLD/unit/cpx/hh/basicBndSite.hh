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

#ifndef CPX_BASICBNDSITE_H
#define CPX_BASICBNDSITE_H

#include <string>
#include <map>
#include <set>
#include "cpx/siteShape.hh"
#include "utl/xcpt.hh"

namespace cpx
{
  class basicBndSite
  {
  protected:
    std::string name;

    const siteShape* pDefaultShape;

  public:
    // Throws an exception if the name for the default shape
    // is not one of the shape names.
    basicBndSite(const std::string& rName,
		 const std::set<std::string>& rShapeNames,
		 const std::string& rDefaultShapeName)
      throw(utl::xcpt);

    // This is a hack that I would like to eliminate: resetting
    // the default shape pointer.
    //
    // The only alternatives I see now are to make the shapes "permanent"
    // (i.e. managed) or binding sites themselves, or both.
    //
    // Fortunately, copying binding sites is done only once or twice.
    basicBndSite(const basicBndSite& rOriginal)
      throw(utl::xcpt);
    
    const std::string&
    getName(void) const
    {
      return name;
    }

    // These have to be traversed in parsing dimerizations; hence
    // exposed publicly.
    std::map<std::string, siteShape> shapesByName;

    // Returns null if site has no shape with the given name.
    const siteShape*
    getShape(const std::string& rShapeName) const;

    // Throws exception if site has no shape with the given name.
    const siteShape*
    mustGetShape(const std::string& rShapeName)
      throw(utl::xcpt);

    const siteShape*
    getDefaultShape(void) const
    {
      return pDefaultShape;
    }
  };
}

#endif // CPX_BASICBNDSITE_H
