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

#ifndef PLEXMAP_H
#define PLEXMAP_H

/*! \file plexMap.hh
  \ingroup plexStructGroup
  \brief Defines maps, isomorphisms of complexes.

  A map from one plex to another is a map from mols to mols and
  bindings to bindings such that corresponding bindings bind
  corresponding mols at the same binding site indices.

  An isomorphism, or equivalence, between plexes is a pair of plex
  maps, each inverse to the other.  When there is an isomorphism
  between two plexes, they are in the same structural class; that is,
  they are the same complex, but possibly with the mols and bindings
  in different positions.*/

#include <vector>
#include <iostream>
#include "plex/plex.hh"

namespace plx
{
  /*! \ingroup plexStructGroup
    \brief Structure-preserving map from one complex to another. */
  class plexMap
  {
  public:
    std::vector<int> molMap;
    std::vector<int> bindingMap;

    plexMap(void)
    {}

    plexMap(int molCount,
	    int bindingCount) :
      molMap(molCount, -1),
      bindingMap(bindingCount, -1)
    {}

    // Methods for applying a plexMap to the specs of various physical
    // features of the complex.
    plexSiteSpec
    applyToSiteSpec(const plexSiteSpec& rSourceSiteSpec) const
    {
      plexSiteSpec targetSpec;
      targetSpec.setMolNdx(molMap[rSourceSiteSpec.molNdx()]);
      targetSpec.setSiteNdx(rSourceSiteSpec.siteNdx());
      return targetSpec;
    }

    // Ever used?
    plexBindingSpec
    applyToBindingSpec(plexBindingSpec sourceBindingSpec) const
    {
      return bindingMap[sourceBindingSpec];
    }

    // Ever used?
    plexMolSpec
    applyToMolSpec(plexMolSpec sourceMolSpec) const
    {
      return molMap[sourceMolSpec];
    }

    // Leaving out subPlexSpec for the time being; it's slightly
    // more complicated, and I doubt it's needed.

    bool molUnmapped(int molNdx)
    {
      return molMap[molNdx] < 0;
    }

    bool bindingUnmapped(int bindingNdx)
    {
      return bindingMap[bindingNdx] < 0;
    }

    bool canMapMol(const plex& rSrcPlex,
		   int srcMolNdx,
		   const plex& rTgtPlex,
		   int tgtMolNdx);

    bool canMapSite(const plex& rSrcPlex,
		    const plexSiteSpec& rSrcSite,
		    const plex& rTgtPlex,
		    const plexSiteSpec& rTgtSite);

    bool canMapBinding(const plex& rSrcPlex,
		       int fromBindingNdx,
		       const plex& rTgtPlex,
		       int toBindingNdx,
		       bool& rMustFlip);
  
    void
    doMapBinding(const plex& rSrcPlex,
		 int srcBindingNdx,
		 const plex& rTgtPlex,
		 int tgtBindingNdx,
		 bool flipBinding);

    static plexMap
    makeIdentity(int molCount,
		 int bindingCount)
    {
      plexMap returnValue(molCount,
			  bindingCount);
      while(0 < molCount--)
	returnValue.molMap[molCount] = molCount;
    
      while(0 < bindingCount--)
	returnValue.bindingMap[bindingCount] = bindingCount;

      return returnValue;
    }

    // So that plexMaps and plexIsoPairs can be sorted/used for map keys.
    bool
    operator<(const plexMap& rRightMap) const
    {
      return ((molMap < rRightMap.molMap)
	      || ((molMap == rRightMap.molMap)
		  && (bindingMap < rRightMap.bindingMap)));
    }
    bool
    operator==(const plexMap& rRightMap) const
    {
      return ((molMap == rRightMap.molMap)
	      && (bindingMap == rRightMap.bindingMap));
    }
  };

  std::ostream& operator<<(std::ostream& rOstr, plexMap& rIso);

  /*! \ingroup plexStructGroup
    \brief An isomorphism of complexes.

    When a plexIsoPair exists between two complexes, then they are the
    same complex, but with the mols and bindings scrambled up.  The
    plexIsoPair gives the "unscrambling" maps, in both directions. */
  class plexIsoPair
  {
  public:
    plexMap forward;
    plexMap backward;

    plexIsoPair(void)
    {}

    plexIsoPair(int molCount,
		int bindingCount) :
      forward(molCount, bindingCount),
      backward(molCount, bindingCount)
    {}

    // This is for where plexIsoPair is used to represent an
    // injection.
    plexIsoPair(int forwardMolCount,
		int forwardBindingCount,
		int backwardMolCount,
		int backwardBindingCount) :
      forward(forwardMolCount,
	      forwardBindingCount),
      backward(backwardMolCount,
	       backwardBindingCount)
    {}
  

    bool
    tryMapBinding(const plex& rSrcPlex,
		  int srcBindingNdx,
		  const plex& rTgtPlex,
		  int tgtBindingNdx);

    static plexIsoPair
    makeIdentity(int molCount,
		 int bindingCount)
    {
      plexIsoPair returnValue;
      returnValue.forward = plexMap::makeIdentity(molCount,
						  bindingCount);
      returnValue.backward = plexMap::makeIdentity(molCount,
						   bindingCount);
      return returnValue;
    }

    // So that plexIsoPairs can be used as map keys (parameters for
    // subPlexes.)
    //
    // It seems like I should be able to get away with not comparing
    // one of the maps, since they are supposed to be inverse to one
    // another.
    bool
    operator<(const plexIsoPair& rRightPair) const
    {
      return ((forward < rRightPair.forward)
	      || ((forward == rRightPair.forward)
		  && backward < rRightPair.backward));
    }
  };

  std::ostream& operator<<(std::ostream& rOstr, plexIsoPair& rIso);

  /*! \ingroup plexStructGroup
    \brief Algorithm to determine when complexes are structurally the same.

    This version does not report the identifying plexIsoPair if it is
    found */
  class plexIsoSearch
  {
    const plex& rLeft;
    const plex& rRight;

    bool
    mapRestBindings(int leftBindingIndex,
		    const plexIsoPair& rCurrentIso) const;
  public:
    plexIsoSearch(const plex& rLeftPlex,
		  const plex& rRightPlex) :
      rLeft(rLeftPlex),
      rRight(rRightPlex)
    {}

    virtual
    ~plexIsoSearch(void)
    {}

    virtual void
    onSuccess(const plexIsoPair& rIso) const
    {}

    bool findIso(void) const;

    bool findInjection(void) const;
  };

  /*! \ingroup plexStructGroup

  \brief Algorithm to determine when complexes are structurally the
  same.

  This version reports the identifying plexIsoPair if it is found. */
  class reportIsoSearch : public plexIsoSearch
  {
    plexIsoPair& rReport;
  public:
    reportIsoSearch(const plex& rLeftPlex,
		    const plex& rRightPlex,
		    plexIsoPair& rReportIso) :
      plexIsoSearch(rLeftPlex,
		    rRightPlex),
      rReport(rReportIso)
    {}

    void
    onSuccess(const plexIsoPair& rIso) const
    {
      rReport = rIso;
    }
  };
}

#endif
