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

#ifndef CPX_PLEXISOIMPL_H
#define CPX_PLEXISOIMPL_H

namespace cpx
{
  template<class plexT>
  bool
  plexIso::
  tryMapBinding(const plexT& rSrcPlex,
		int srcBindingNdx,
		const plexT& rTgtPlex,
		int tgtBindingNdx)
  {
    bool flipForward, flipBackward;

    // Test both the forward and backward maps for compatibility with
    // the given mapping of bindings, noting whether the binding needs
    // to be flipped or not.
    //
    // If the binding needs to be flipped for the forward map, then it
    // should need to be flipped for the backward map, too.  This shows
    // that I'm doing the flip check twice, but it's really a cost-free
    // byproduct the necessary check of consistency of the mol mapping,
    // I think...
    if(forward.canMapBinding(rSrcPlex,
			     srcBindingNdx,
			     rTgtPlex,
			     tgtBindingNdx,
			     flipForward)
       &&
       backward.canMapBinding(rTgtPlex,
			      tgtBindingNdx,
			      rSrcPlex,
			      srcBindingNdx,
			      flipBackward))
      {
	forward.doMapBinding(rSrcPlex,
			     srcBindingNdx,
			     rTgtPlex,
			     tgtBindingNdx,
			     flipForward);

	backward.doMapBinding(rTgtPlex,
			      tgtBindingNdx,
			      rSrcPlex,
			      srcBindingNdx,
			      flipBackward);
	return true;
      }
    else return false;
  }
}

#endif // CPX_PLEXISOIMPL_H
