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

#include "mzr/util.hh"
#include "mol/siteShape.hh"
#include "plex/plexSpec.hh"
#include "plex/plexMap.hh"
#include "plex/alloSiteQuery.hh"
#include "plex/plexQuery.hh"

namespace plx
{
  bnd::siteParam
  siteToShapeMap::getSiteShape(const plexSiteSpec& rPlexSiteSpec) const
    throw(unmappedSiteSpecXcpt)
  {
    std::map<plexSiteSpec, bnd::siteParam>::const_iterator iEntry
      = find(rPlexSiteSpec);

    if(iEntry == end())
      {
	throw unmappedSiteSpecXcpt();
      }
    else
      {
	return iEntry->second;
      }
  }

  void
  siteToShapeMap::setSiteShape(const plexSiteSpec& rPlexSiteSpec,
			       const bnd::siteParam& rSiteParam)
  {
    mzr::forceInsert(*this,
		     std::pair<plexSiteSpec, bnd::siteParam>(rPlexSiteSpec,
							     rSiteParam));
  }

  void
  siteToShapeMap::setSiteShape(const plexSiteSpec& rPlexSiteSpec,
			       const bnd::siteParam& rSiteParam,
			       const subPlexSpec& rSubPlexSpec)
  {
    const plexIsoPair& rInjection = rSubPlexSpec.getInjection();
    plexSiteSpec mappedSpec
      = rInjection.forward.applyToSiteSpec(rPlexSiteSpec);
    
    setSiteShape(mappedSpec,
		 rSiteParam);
  }

  void
  siteToShapeMap::setSiteShapes(const siteToShapeMap& rSiteToShapeMap)
  {
    mzr::forceInsert(*this,
		     rSiteToShapeMap.begin(),
		     rSiteToShapeMap.end());
  }

  // Auxiliary function for siteToShapeMap::setSiteShapes below.
  class setSiteShapeTracked :
    public std::unary_function<siteToShapeMap::value_type, void>
  {
    siteToShapeMap& rTarget;
    const subPlexSpec& rSpec;
  public:
    setSiteShapeTracked(siteToShapeMap& rTargetMap,
			const subPlexSpec& rSubPlexSpec) :
      rTarget(rTargetMap),
      rSpec(rSubPlexSpec)
    {}

    void
    operator()(const argument_type& rSiteShapePair) const
    {
      rTarget.setSiteShape(rSiteShapePair.first,
			   rSiteShapePair.second,
			   rSpec);
    }
  };

  void
  siteToShapeMap::setSiteShapes(const siteToShapeMap& rSiteToShapeMap,
				const subPlexSpec& rSubPlexSpec)
  {
    std::for_each(rSiteToShapeMap.begin(),
		  rSiteToShapeMap.end(),
		  setSiteShapeTracked(*this,
				      rSubPlexSpec));
  }

  void
  queryAllosteryList::
  addQueryAndMap(const andPlexQueries* pQuery,
		 const siteToShapeMap& rSiteToShapeMap)
  {
    push_back
      (std::pair<const andPlexQueries*, siteToShapeMap>(pQuery,
							rSiteToShapeMap));
  }

  // Auxiliary function for queryAllosteryList::setSatisfiedQuerySiteShapes
  // below.
  class setShapesIfSatisfied :
    public std::unary_function<queryAllosteryList::value_type, void>
  {
    plexParam& rParam;
    
  public:
    setShapesIfSatisfied(plexParam& rPlexParam) :
      rParam(rPlexParam)
    {}

    void
    operator()(const queryAllosteryList::value_type& rQueryShapeMapPair) const
    {
      const andPlexQueries& rQuery = *(rQueryShapeMapPair.first);
      if(rQuery(rParam))
      {
	const siteToShapeMap& rSiteToShapeMap = rQueryShapeMapPair.second;
	rParam.siteParams.setSiteShapes(rSiteToShapeMap);
      }
    }
  };

  void
  queryAllosteryList::setSatisfiedQuerySiteShapes(plexParam& rPlexParam) const
  {
    std::for_each(begin(),
		  end(),
		  setShapesIfSatisfied(rPlexParam));
  }

  // Auxiliary function for injectSatisfiedQuerySiteShapes below.
  class setShapesIfSatisfiedTracked :
    public std::unary_function<queryAllosteryList::value_type, void>
  {
    plexParam& rParam;
    const subPlexSpec& rSpec;
  public:
    setShapesIfSatisfiedTracked(plexParam& rPlexParam,
				const subPlexSpec& rSubPlexSpec) :
      rParam(rPlexParam),
      rSpec(rSubPlexSpec)
    {}

    void
    operator()(const argument_type& rQueryShapeMapPair) const
    {
      const andPlexQueries& rQuery = *(rQueryShapeMapPair.first);
      if(rQuery.applyTracked(rParam,
			     rSpec))
	{
	  const siteToShapeMap& rSiteToShapeMap
	    = rQueryShapeMapPair.second;
	  rParam.siteParams.setSiteShapes(rSiteToShapeMap,
					  rSpec);
	}
    }
  };

  void
  queryAllosteryList::
  setSatisfiedQuerySiteShapes(plexParam& rPlexParam,
			      const subPlexSpec& rSubPlexSpec) const
  {
    std::for_each(begin(),
		  end(),
		  setShapesIfSatisfiedTracked(rPlexParam,
					      rSubPlexSpec));
  }
}
