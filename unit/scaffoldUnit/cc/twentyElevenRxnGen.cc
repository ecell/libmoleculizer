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

#include "plex/cxOmniParam.hh"
#include "scaffold/twentyElevenRxnGen.hh"

namespace scaf
{
  void
  twentyElevenRxnGen::makeReactions(const plx::omniInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the notifying complex, which contains the "enabling"
	// subcomplex, in accessors.
	plx::cxOmni cxEnabler(rContext);

	// Get the state of Ste20, to determine if Ste20 is bound to ATP.
	//
	// First, translate the index of Ste20 in the enabling omniplex
	// into its index in the notifying complex, using the injection
	// of the omniplex into the notifying complex that was created
	// when the notifying complex was recognized as containing the
	// enabling subcomplex.
	plx::plexMolSpec ste20SpecTr = cxEnabler.translateMolSpec(ste20Spec);
	// Get the state of Ste20 in the notifying complex.
	bnd::molParam ste20Prm = cxEnabler.getMolParam(ste20SpecTr);
	const bnd::modMolState& rSte20State = pSte20->externState(ste20Prm);

	// Generate new reactions only if Ste20 is bound to ATP.
	if(rSte20State[atpModSiteNdx] == pAtpBound)
	  {
	    // Get the state of Ste11, in order to see which of the
	    // target phosphorylation sites are open.
	    plx::plexMolSpec ste11SpecTr
	      = cxEnabler.translateMolSpec(ste11Spec);
	    bnd::molParam ste11Prm = cxEnabler.getMolParam(ste11SpecTr);
	    const bnd::modMolState& rSte11State
	      = pSte11->externState(ste11Prm);

	    // Generate a phosphorylation reaction for each of the modification
	    // sites indicated by the "activity mask" if the site is free.
	    //
	    // We need indices here.  The "official" way would be to
	    // use a generator in conjunction with the non-existent
	    // binary version of std::for_each.
	    int targetModSiteNdx = activityMask.size();
	    while(0 < targetModSiteNdx--)
	      {
		// Is the modification site "active", i.e. specified as a
		// phosphorylation target of Ste20, and is it free?
		if(activityMask[targetModSiteNdx]
		   && (rSte11State[targetModSiteNdx] == pNone))
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

		    // Make the new state of Ste11 obtained by phosphorylating
		    // the free modification site.
		    bnd::modMolState nuSte11State(rSte11State);
		    nuSte11State[targetModSiteNdx] = pPhosphorylated;
		    bnd::molParam nuSte11Prm
		      = pSte11->internState(nuSte11State);

		    // Make the new state of Ste20 by changing the modification
		    // from atp-bound to adp-bound.
		    bnd::modMolState nuSte20State(rSte20State);
		    nuSte20State[atpModSiteNdx] = pAdpBound;
		    bnd::molParam nuSte20Prm
		      = pSte20->internState(nuSte20State);

		    // Assemble the new molparams of the notifying complex.
		    std::vector<bnd::molParam>
		      nuMolParams(cxEnabler.getMolParams());
		    nuMolParams[ste11SpecTr] = nuSte11Prm;
		    nuMolParams[ste20SpecTr] = nuSte20Prm;

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
}
