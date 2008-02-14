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

#include <sstream>
#include "cpx/basicBndSite.hh"

namespace cpx
{
  class unknownSiteShapeXcpt : 
    public utl::xcpt
  {
    static std::string
    mkMsg(const basicBndSite& rSite,
	  const std::string& rBadSiteShapeName)
    {
      std::ostringstream msgStream;
      msgStream << "Binding site "
		<< rSite.getName()
		<< " has no shape named "
		<< rBadSiteShapeName
		<< ".";
      return msgStream.str();
    }
  public:
    unknownSiteShapeXcpt(const basicBndSite& rSite,
			 const std::string& rBadSiteShapeName) :
      utl::xcpt(mkMsg(rSite,
		      rBadSiteShapeName))
    {}
  };

  class insertShape :
    public std::unary_function<std::string, void>
  {
    std::map<std::string, siteShape>& rShapeMap;

  public:
    insertShape(std::map<std::string, siteShape>& rShapesByName) :
      rShapeMap(rShapesByName)
    {}

    void
    operator()(const std::string& rShapeName)
    {
      rShapeMap.insert
	(std::pair<std::string, siteShape>(rShapeName,
					   siteShape(rShapeName)));
    }
  };

  basicBndSite::
  basicBndSite(const std::string& rName,
	       const std::set<std::string>& rShapeNames,
	       const std::string& rDefaultShapeName) 
    throw(utl::xcpt) :
    name(rName)
  {
    std::for_each(rShapeNames.begin(),
		  rShapeNames.end(),
		  insertShape(shapesByName));

    pDefaultShape = mustGetShape(rDefaultShapeName);
  }

  basicBndSite::
  basicBndSite(const basicBndSite& rOriginal)
    throw(utl::xcpt) :
    name(rOriginal.getName()),
    shapesByName(rOriginal.shapesByName)
  {
    pDefaultShape = mustGetShape(rOriginal.getDefaultShape()->getName());
  }

  const siteShape*
  basicBndSite::
  getShape(const std::string& rShapeName) const
  {
    std::map<std::string, siteShape>::const_iterator iShapeEntry
      = shapesByName.find(rShapeName);

    return (shapesByName.end() == iShapeEntry)
      ? 0
      : &(iShapeEntry->second);
  }

  const siteShape*
  basicBndSite::
  mustGetShape(const std::string& rShapeName)
    throw(utl::xcpt)
  {
    const siteShape* pShape
      = getShape(rShapeName);

    if(! pShape) throw unknownSiteShapeXcpt(*this,
					    rShapeName);
    return pShape;
  }
}


