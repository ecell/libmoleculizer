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

#include "cpx/unmappedSiteSpecXcpt.hh"

namespace cpx
{
  template<class omniPlexT>
  void
  siteToShapeMap::
  setSiteShape(const siteSpec& rPlexSiteSpec,
	       const siteShape* pSiteShape,
	       const subPlexSpec<omniPlexT>& rSubPlexSpec)
  {
    const plexIso& rInjection
      = rSubPlexSpec.getInjection();

    siteSpec mappedSpec
      = rInjection.forward.applyToSiteSpec(rPlexSiteSpec);
    
    setSiteShape(mappedSpec,
		 pSiteShape);
  }

  template<class omniPlexT>
  class setSiteShapeTracked :
    public std::unary_function<siteToShapeMap::value_type, void>
  {
  public:
    typedef subPlexSpec<omniPlexT> subPlexSpecType;
    
    siteToShapeMap& rTarget;
    const subPlexSpecType& rSpec;
  public:
    setSiteShapeTracked(siteToShapeMap& rTargetMap,
			const subPlexSpecType& rSubPlexSpec) :
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

  template<class omniPlexT>
  void
  siteToShapeMap::
  setSiteShapes(const siteToShapeMap& rSiteToShapeMap,
		const subPlexSpec<omniPlexT>& rSubPlexSpec)
  {
    std::for_each(rSiteToShapeMap.begin(),
		  rSiteToShapeMap.end(),
		  setSiteShapeTracked<omniPlexT>(*this,
						 rSubPlexSpec));
  }
}
