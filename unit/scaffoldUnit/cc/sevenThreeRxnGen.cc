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
#include "scaffold/sevenThreeRxnGen.hh"

namespace scaf
{
  void
  sevenThreeRxnGen::makeReactions(const plx::omniInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the notifying complex, which contains the enabling omniplex,
	// in accessors.
	plx::cxOmni cxEnabler(rContext);

	// Get the state of Ste7, in order to determine if it is bound to ATP
	// and if it's phosphorylation state permits it to be active.
	plx::plexMolSpec ste7SpecTr = cxEnabler.translateMolSpec(ste7Spec);
	bnd::molParam ste7Prm = cxEnabler.getMolParam(ste7SpecTr);
	const bnd::modMolState& rSte7State = pSte7->externState(ste7Prm);

	// Generate new reaction only if Ste7 is bound to ATP, is
	// phosphorylated at the "activating" phosphorylation site targeted by
	// Ste11, and is NOT phosphorylated at the "inhibiting"
	// phosphorylation site targeted by Fus3.
	if((rSte7State[atpModNdx] == pAtpBound)
	   && (rSte7State[activeModNdx] == pPhosphorylated)
	   && (rSte7State[inhibModNdx] == pNone))
	  {
	    // Get the state of Fus3, in order to determine which
	    // phosphorylation sites available to Ste7 are free.
	    plx::plexMolSpec fus3SpecTr = cxEnabler.translateMolSpec(fus3Spec);
	    bnd::molParam fus3Prm = cxEnabler.getMolParam(fus3SpecTr);
	    const bnd::modMolState& rFus3State = pFus3->externState(fus3Prm);

	    // Create phosphorylation reaction for each available target site.
	    // This requires the use of indices.
	    int modSiteNdx = activityMask.size();
	    while(0 < modSiteNdx--)
	      {
		// Is the modification site available for phosphorylation
		// and is it free?
		if(activityMask[modSiteNdx]
		   && (rFus3State[modSiteNdx] == pNone))
		  {
		    // Create phosphorylation reaction.
		    mzr::reaction* pReaction = new mzr::reaction();
		    pFamily->addReaction(pReaction,
					 rMzrUnit);

		    // Add the notifying complex as the only substrate of the
		    // reaction.
		    pReaction->addSubstrate(cxEnabler.getSpecies(),
					    1);
		    pReaction->sensitizeToSubstrates(rMzrUnit);

		    // Extrapolate the reaction rate.
		    pReaction->setRate(pExtrap->getRate(rContext));

		    // Make the new state of Fus3 obtained by phosphorylating
		    // the free modification site.
		    bnd::modMolState nuFus3State(rFus3State);
		    nuFus3State[modSiteNdx] = pPhosphorylated;
		    bnd::molParam nuFus3Prm = pFus3->internState(nuFus3State);

		    // Make the new state of Ste7 obtained by changing the
		    // modification from atp-bound to adp-bound.
		    bnd::modMolState nuSte7State(rSte7State);
		    nuSte7State[atpModNdx] = pAdpBound;
		    bnd::molParam nuSte7Prm = pSte7->internState(nuSte7State);

		    // Assemble the new molparams of the notifying complex.
		    std::vector<bnd::molParam>
		      nuMolParams(cxEnabler.getMolParams());
		    nuMolParams[fus3SpecTr] = nuFus3Prm;
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
}

