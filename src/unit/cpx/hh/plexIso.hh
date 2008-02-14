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

#ifndef CPX_PLEXISO_H
#define CPX_PLEXISO_H

#include "cpx/plexMap.hh"

namespace cpx
{
  class plexIso
  {
  public:
    plexMap forward;
    plexMap backward;

    plexIso(void)
    {}

    plexIso(int molCount,
		int bindingCount) :
      forward(molCount, bindingCount),
      backward(molCount, bindingCount)
    {}

    // This is for where plexIso is used to represent an
    // injection.
    plexIso(int forwardMolCount,
		int forwardBindingCount,
		int backwardMolCount,
		int backwardBindingCount) :
      forward(forwardMolCount,
	      forwardBindingCount),
      backward(backwardMolCount,
	       backwardBindingCount)
    {}
  

    // Tests if the isomorphism can be extended by the given assignment of
    // bindings, either with or without the "flip." If the extension is
    // possible, it is performed on 'this' isomorphism, and true is returned.
    // Otherwise, false is returned.
    template<class plexT>
    bool
    tryMapBinding(const plexT& rSrcPlex,
		  int srcBindingNdx,
		  const plexT& rTgtPlex,
		  int tgtBindingNdx);


    // Generates the identity isomorphism for a plex with the given
    // number of mols and bindings.
    static plexIso
    makeIdentity(int molCount,
		 int bindingCount);

    // So that plexIsos can be used as map keys (parameters for
    // subPlexes.)
    //
    // It seems like I should be able to get away with not comparing
    // one of the maps, since they are supposed to be inverse to one
    // another.
    bool
    operator<(const plexIso& rRightIso) const
    {
      return ((forward < rRightIso.forward)
	      || ((forward == rRightIso.forward)
		  && backward < rRightIso.backward));
    }
  };
}

#include "cpx/plexIsoImpl.hh"

#endif // CPX_PLEXISO_H
