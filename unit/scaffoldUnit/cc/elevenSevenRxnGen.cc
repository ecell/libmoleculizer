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

#include "mzr/mzrUnit.hh"
#include "plex/cxOmniParam.hh"
#include "scaffold/elevenSevenRxnGen.hh"

namespace scaf
{
  void
  elevenSevenRxnGen::makeReactions(const plx::omniInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the notifying complex, which contains the enabling complex,
	// in accessors.
	plx::cxOmni cxEnabler(rContext);

	// Get the state of Ste11, in order to determine if it is bound to ATP
	// and is phosphorylated the requisite number of times.
	plx::plexMolSpec ste11SpecTr = cxEnabler.translateMolSpec(ste11Spec);
	bnd::molParam ste11Prm = cxEnabler.getMolParam(ste11SpecTr);
	const bnd::modMolState& rSte11State = pSte11->externState(ste11Prm);

	// Generate new reaction only if Ste11 is bound to ATP and
	// phosphorylated sufficiently many times.
	if((rSte11State[atpModSiteNdx] == pAtpBound)
	   && (elevenActiveCount <= rSte11State.modCount(pPhosphorylated)))
	  {
	    // Get the state of Ste7, in order to determine if Ste11's target
	    // phosphorylation site is free.
	    plx::plexMolSpec ste7SpecTr = cxEnabler.translateMolSpec(ste7Spec);
	    bnd::molParam ste7Prm = cxEnabler.getMolParam(ste7SpecTr);
	    const bnd::modMolState& rSte7State = pSte7->externState(ste7Prm);

	    // Generate new reaction only if phosphorylation site is free.
	    if(rSte7State[targetModSiteNdx] == pNone)
	      {
		// Construct a phosphorylation reaction.
		mzr::reaction* pReaction = new mzr::reaction();
		pFamily->addReaction(pReaction,
				     rMzrUnit);

		// Add the complex as the only substrate of the reaction.
		pReaction->addSubstrate(cxEnabler.getSpecies(),
					1);
		pReaction->sensitizeToSubstrates(rMzrUnit);

		// Set the (unary) reaction rate.
		pReaction->setRate(pExtrap->getRate(rContext));

		// Make the new state of Ste11 by changing the modification
		// from atp-bound to adp-bound.
		bnd::modMolState nuSte11State(rSte11State);
		nuSte11State[atpModSiteNdx] = pAdpBound;
		bnd::molParam nuSte11Prm = pSte11->internState(nuSte11State);

		// Make the new state of Ste7 obtained by phosphorylating
		// the target modification site.
		bnd::modMolState nuSte7State(rSte7State);
		nuSte7State[targetModSiteNdx] = pPhosphorylated;
		bnd::molParam nuSte7Prm = pSte7->internState(nuSte7State);

		// Assemble the new molparams of the notifying complex.
		std::vector<bnd::molParam> nuMolParams(cxEnabler.getMolParams());
		nuMolParams[ste11SpecTr] = nuSte11Prm;
		nuMolParams[ste7SpecTr] = nuSte7Prm;

		// Make the product species.
		plx::plexSpecies* pResult
		  = cxEnabler.getPlexFamily().getMember(nuMolParams);

		// Add it as a product of multiplicity 1.
		pReaction->addProduct(pResult,
				      1);
	      }
	  }
      }
  }
}
