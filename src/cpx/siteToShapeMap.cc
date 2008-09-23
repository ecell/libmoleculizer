/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008 The Molecular Sciences Institute.
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001

//                     Email: lok@molsci.org
//   
/////////////////////////////////////////////////////////////////////////////

#include "utl/forceInsert.hh"
#include "cpx/siteToShapeMap.hh"
#include "cpx/subPlexSpec.hh"
#include "cpx/unmappedSiteSpecXcpt.hh"

namespace cpx
{
  const siteShape*
  siteToShapeMap::
  getSiteShape(const siteSpec& rPlexSiteSpec) const
  {
    const_iterator iEntry 
      = find(rPlexSiteSpec);

    return end() == iEntry
      ? 0
      : iEntry->second;
  }
  
  const siteShape*
  siteToShapeMap::
  mustGetSiteShape(const siteSpec& rPlexSiteSpec) const
    throw(utl::xcpt)
  {
    const siteShape* pSiteShape
      = getSiteShape(rPlexSiteSpec);
    
    if(! pSiteShape) throw unmappedSiteSpecXcpt();

    return pSiteShape;
  }

  void
  siteToShapeMap::
  setSiteShape(const siteSpec& rPlexSiteSpec,
	       const siteShape* pSiteShape)
  {
    utl::forceInsert(*this,
		     std::pair<siteSpec, const siteShape*>(rPlexSiteSpec,
							   pSiteShape));
  }

  class setOneSiteShape :
    public std::unary_function<siteToShapeMap::value_type, void>
  {
    siteToShapeMap& rTargetMap;
  public:
    setOneSiteShape(siteToShapeMap& refTargetMap) :
      rTargetMap(refTargetMap)
    {}

    void
    operator()(const argument_type& rSpecShapePair) const
    {
      rTargetMap.setSiteShape(rSpecShapePair.first,
			      rSpecShapePair.second);
    }
  };

  void
  siteToShapeMap::
  setSiteShapes(const siteToShapeMap& rSiteToShapeMap)
  {
    std::for_each(rSiteToShapeMap.begin(),
		  rSiteToShapeMap.end(),
		  setOneSiteShape(*this));
  }
}
