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

#include "plex/plex.hh"
#include "plex/plexSpec.hh"
#include "plex/plexMap.hh"
#include "mzr/moleculizer.hh"

namespace plx
{
  bool
  plexMap::canMapMol(const plex& rSrcPlex,
		     int srcMolNdx,
		     const plex& rTgtPlex,
		     int tgtMolNdx)
  {
    return 
      (molMap[srcMolNdx] == tgtMolNdx)
      || (molUnmapped(srcMolNdx) && (rSrcPlex.mols[srcMolNdx]
				     == rTgtPlex.mols[tgtMolNdx]));
  }

  bool
  plexMap::canMapSite(const plex& rSrcPlex,
		      const plexSiteSpec& rSrcSite,
		      const plex& rTgtPlex,
		      const plexSiteSpec& rTgtSite)
  {
    return
      canMapMol(rSrcPlex,
		rSrcSite.molNdx(),
		rTgtPlex,
		rTgtSite.molNdx())
      &&
      rSrcSite.siteNdx() == rTgtSite.siteNdx();
  }

  //  bool
  //  plexMap::canMapBinding(const plex& rSrcPlex,
  //  		       int fromBindingNdx,
  //  		       const plex& rTgtPlex,
  //  		       int toBindingNdx,
  //  		       bool& rMustFlip)
  //  {
  //    const plexSiteSpec& rLeftSrcSite
  //      = rSrcPlex.bindings[fromBindingNdx].leftSite();
  //    const plexSiteSpec& rRightSrcSite
  //      = rSrcPlex.bindings[fromBindingNdx].rightSite();
  //    const plexSiteSpec& rLeftTgtSite
  //      = rTgtPlex.bindings[toBindingNdx].leftSite();
  //    const plexSiteSpec& rRightTgtSite
  //      = rTgtPlex.bindings[toBindingNdx].rightSite();

  //    if(canMapSite(rSrcPlex,
  //  		rLeftSrcSite,
  //  		rTgtPlex,
  //  		rLeftTgtSite)
  //       &&
  //       canMapSite(rSrcPlex,
  //  		rRightSrcSite,
  //  		rTgtPlex,
  //  		rRightTgtSite))
  //      {
  //        rMustFlip = false;
  //      }
  //    else if(canMapSite(rSrcPlex,
  //  		     rRightSrcSite,
  //  		     rTgtPlex,
  //  		     rLeftTgtSite)
  //  	  &&
  //  	  canMapSite(rSrcPlex,
  //  		     rLeftSrcSite,
  //  		     rTgtPlex,
  //  		     rRightTgtSite))
  //      {
  //        rMustFlip = true;
  //      }
  //    else return false;

  //    // I can't see at this time (Thu Jul 26 11:44:49 PDT 2001) why I
  //    // didn't make this a precondition to the above checking of the site
  //    // map.  This would also have the advantage of not setting rMustFlip
  //    // unless a mapping is found.  Could I have been using the rMustFlip
  //    // value while ignoring the mapping status of the binding somewhere?
  //    // In other words, could this have been intentional?
  //    return bindingUnmapped(fromBindingNdx)
  //      || bindingMap[fromBindingNdx] == toBindingNdx;
  //  }

  // This routine returns true and gives (modifies) the flag rMustFlip
  // if the given mapping of bindings is consistent with the current
  // state of the plexMap.  The rMustFlip flag is determined even if the
  // given mapping of bindings is consistent because it is already in
  // the map.  If the given mapping is not consistent with the current
  // state of the plexMap, then rMustFlip is left unmodified and false
  // is returned.
  bool
  plexMap::canMapBinding(const plex& rSrcPlex,
			 int fromBindingNdx,
			 const plex& rTgtPlex,
			 int toBindingNdx,
			 bool& rMustFlip)
  {
    // First check the binding map to see if there is already
    // a conflicting mapping of the binding.
    if(bindingUnmapped(fromBindingNdx)
       || bindingMap[fromBindingNdx] == toBindingNdx)
      {
	// Since the binding can be mapped, or is already mapped,
	// we check whether the mols can be mapped consistently,
	// noting whether the binding has to be "flipped" to get
	// the consistent mapping of mols.

	const plexSiteSpec& rLeftSrcSite
	  = rSrcPlex.bindings[fromBindingNdx].leftSite();
	const plexSiteSpec& rRightSrcSite
	  = rSrcPlex.bindings[fromBindingNdx].rightSite();
	const plexSiteSpec& rLeftTgtSite
	  = rTgtPlex.bindings[toBindingNdx].leftSite();
	const plexSiteSpec& rRightTgtSite
	  = rTgtPlex.bindings[toBindingNdx].rightSite();

	// Try to match up the mols in the unflipped order.
	if(canMapSite(rSrcPlex,
		      rLeftSrcSite,
		      rTgtPlex,
		      rLeftTgtSite)
	   &&
	   canMapSite(rSrcPlex,
		      rRightSrcSite,
		      rTgtPlex,
		      rRightTgtSite))
	  {
	    // Unflipped mapping of the mols is consistent.
	    rMustFlip = false;
	    return true;
	  }
	// Try to match up the mols in the flipped order.
	else if(canMapSite(rSrcPlex,
			   rRightSrcSite,
			   rTgtPlex,
			   rLeftTgtSite)
		&&
		canMapSite(rSrcPlex,
			   rLeftSrcSite,
			   rTgtPlex,
			   rRightTgtSite))
	  {
	    // Flipped mapping of the mols is consistent.
	    rMustFlip = true;
	    return true;
	  }
	// Although the binding map is okay, there is no consistent
	// way to map the mols.
	else return false;
      }
    // The binding map conflicts with the proposed mapping of the
    // binding.
    else return false;
  }

  // Given that a mapping of a binding is possible, maybe by flipping
  // the mols, this routine actually modifies the map by adding the
  // mapping of the binding and the appropriate flipped or unflipped
  // mapping of the mols.
  void
  plexMap::doMapBinding(const plex& rSrcPlex,
			int srcBindingNdx,
			const plex& rTgtPlex,
			int tgtBindingNdx,
			bool flipBinding)
  {
    int leftSrcMolNdx
      = rSrcPlex.bindings[srcBindingNdx].leftSite().molNdx();
    int rightSrcMolNdx
      = rSrcPlex.bindings[srcBindingNdx].rightSite().molNdx();
    int leftTgtMolNdx
      = rTgtPlex.bindings[tgtBindingNdx].leftSite().molNdx();
    int rightTgtMolNdx
      = rTgtPlex.bindings[tgtBindingNdx].rightSite().molNdx();

    if(flipBinding)
      {
	molMap[leftSrcMolNdx] = rightTgtMolNdx;
	molMap[rightSrcMolNdx] = leftTgtMolNdx;
      }
    else
      {
	molMap[leftSrcMolNdx] = leftTgtMolNdx;
	molMap[rightSrcMolNdx] = rightTgtMolNdx;
      }

    bindingMap[srcBindingNdx] = tgtBindingNdx;
  }

  // If the given map of bindings can be done both forward and backward,
  // this routine does it and returns true.  Otherwise, it does nothing
  // and returns false.
  bool
  plexIsoPair::tryMapBinding(const plex& rSrcPlex,
			     int srcBindingNdx,
			     const plex& rTgtPlex,
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
	// Debugging code.
	if((flipForward && (! flipBackward))
	   ||
	   (flipBackward && (! flipForward)))
	  throw incorrectFlipFlagsXcpt();
      
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

  std::ostream&
  operator<<(std::ostream& rOstr, plexMap& rMap)
  {
    rOstr << "molMap: ";
    for(int molNdx = 0;
	molNdx < (int) rMap.molMap.size();
	molNdx++)
      {
	rOstr << rMap.molMap[molNdx]
	      << ' ';
      }
    rOstr << std::endl
	  << "bindingMap: ";
    for(int bindingNdx = 0;
	bindingNdx < (int) rMap.bindingMap.size();
	bindingNdx++)
      {
	rOstr << rMap.bindingMap[bindingNdx]
	      << ' ';
      }
    return rOstr;
  }

  std::ostream&
  operator<<(std::ostream& rOstr, plexIsoPair& rIso)
  {
    rOstr << "forward:"
	  << std::endl
	  << rIso.forward
	  << std::endl
	  << "backward:"
	  << std::endl
	  << rIso.backward;
    return rOstr;
  }

  bool
  plexIsoSearch::mapRestBindings(int leftBindingNdx,
				 const plexIsoPair& rCurrentIso) const
  {
    // Are we done?  
    if(((int) rLeft.bindings.size()) <= leftBindingNdx)
      {
	onSuccess(rCurrentIso);
	return true;
      }

    // Try to extend the given isomorphism over the given left binding
    // until one is found that can be extended to a full isomorpism.
    for(int rightBindingNdx = 0;
	rightBindingNdx < (int) rRight.bindings.size();
	rightBindingNdx++)
      {
	// Construct a new temporary isomorphism for this trial.
	plexIsoPair tmpIso(rCurrentIso);

	if(tmpIso.tryMapBinding(rLeft,
				leftBindingNdx,
				rRight,
				rightBindingNdx)
	   && mapRestBindings(leftBindingNdx + 1,
			      tmpIso))
	  {
	    return true;
	  }
      }
    return false;
  }

  bool
  plexIsoSearch::findInjection(void) const
  {
    plexIsoPair tmpIso(rLeft.mols.size(),
		       rLeft.bindings.size(),
		       rRight.mols.size(),
		       rRight.bindings.size());

    // Since we're assuming that the left plex is connected,
    // we've got two cases, either every mol is on a binding
    // or there is only one mol.
    if(rLeft.bindings.size() > 0)
      {
	// Search for mappings of all the bindings in the pattern.
	if(mapRestBindings(0,
			   tmpIso))
	  {
	    return true;
	  }
      }
    else if(1 == rLeft.mols.size())
      {
	// Attempt to match with each of the mols in the target complex.
	for(int tgtMolNdx = 0;
	    tgtMolNdx < (int) rRight.mols.size();
	    tgtMolNdx++)
	  {
	    if(tmpIso.forward.canMapMol(rLeft,
					0,
					rRight,
					tgtMolNdx))
	      {
		tmpIso.forward.molMap[0] = tgtMolNdx;
		tmpIso.backward.molMap[tgtMolNdx] = 0;
		onSuccess(tmpIso);
		return true;
	      }
	  }
      }
    else
      {
	throw plexNotConnectedXcpt();
      }
  
    return false;
  }

  // To find and isomorophism, find an injection and apply the
  // pigeonhole principle.
  bool
  plexIsoSearch::findIso(void) const
  {
    if((rLeft.bindings.size() != rRight.bindings.size())
       ||
       (rLeft.mols.size() != rRight.mols.size()))
      {
	return false;
      }
    else
      {
	return findInjection();
      }
  }
}
