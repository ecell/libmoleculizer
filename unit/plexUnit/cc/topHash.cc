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

#include <set>
#include "mzr/linearHash.hh"
#include "mol/mol.hh"
#include "plex/plex.hh"

namespace plx
{
  // Function to hash plexes by "topological sorting."
  class hashMolRec
  {
    const plex& rPlex;
    const std::map<plexSiteSpec, int>& rSiteToBindings;
    std::set<int>& rMolsSeen;
  
  public:
    hashMolRec(const plex& refPlex,
	       const std::map<plexSiteSpec, int>& refSiteToBindings,
	       std::set<int>& refMolsSeen) :
      rPlex(refPlex),
      rSiteToBindings(refSiteToBindings),
      rMolsSeen(refMolsSeen)
    {}
  
    size_t operator()(int molNdx,
		      int depth) const
    {
      mzr::linearHash lh;
      size_t hashValue = 0;

      // Have we seen this mol?
      //
      // Not by pointer, since who knows when we might encounter
      // a particular mol species in a random traversal.
      if(rMolsSeen.find(molNdx) == rMolsSeen.end())
	{
	  // We have now seen the mol that this path goes to.
	  rMolsSeen.insert(molNdx);

	  // Initialize the hash value using the mol pointer and
	  // the depth.
  	  bnd::mol* pMol = rPlex.mols[molNdx];
  	  hashValue = lh(lh((size_t) pMol)
			 + (size_t) depth);

	  // Traverse the sites on this mol.  We can use the ordering
	  // of the sites on the mol.
	  for(int siteNdx = 0;
	      siteNdx < (int) pMol->getSiteCount();
	      siteNdx++)
	    {
	      plexSiteSpec siteSpec(molNdx, siteNdx);

	      std::map<plexSiteSpec, int>::const_iterator iPair
		= rSiteToBindings.find(siteSpec);

	      if(iPair != rSiteToBindings.end())
		{
		  int bindingNdx = iPair->second;
		  const plexBinding& rBinding
		    = rPlex.bindings[bindingNdx];

		  if(rBinding.leftSite().molNdx() == molNdx
		     && rBinding.leftSite().siteNdx() == siteNdx)
		    {
		      hashValue
			= lh((size_t) rBinding.rightSite().siteNdx()
			     + lh((size_t) siteNdx
				  + lh(hashValue)));

		      hashValue = lh((*this)(rBinding.rightSite().molNdx(),
					     depth + 1)
				     + lh(hashValue));
		    }
		  else
		    {
		      hashValue
			= lh((size_t) rBinding.leftSite().siteNdx()
			     + lh((size_t) siteNdx
				  + lh(hashValue)));

		      hashValue = lh((*this)(rBinding.leftSite().molNdx(),
					     depth + 1)
				     + lh(hashValue));
		    }

		}
	    }
	}
      return hashValue;
    }
  };

  int
  plex::hashValue(void) const
  {
    std::map<plexSiteSpec, int> siteToBindings;
    makeSiteToBindings(siteToBindings);

    int molNdx = mols.size();
    size_t hashValue = 0;
    while(0 < molNdx--)
      {
	std::set<int> molsSeen;

	hashMolRec hmr(*this,
		       siteToBindings,
		       molsSeen);

	hashValue += hmr(molNdx,
			 0);
      }
    return hashValue;
  }
}
