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
#include "scaffold/threeSevenRxnGen.hh"

namespace scaf
{
  void
  threeSevenRxnGen::makeReactions(const plx::omniInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the notifying complex, which must contain the enabling
	// omniplex, in accessors.
	plx::cxOmni cxEnabler(rContext);

	// Get the state of Fus3, in order to determine if it is bound to ATP
	// and is not phosphorylated more than the limiting number of times.
	plx::plexMolSpec fus3SpecTr = cxEnabler.translateMolSpec(fus3Spec);
	bnd::molParam fus3Prm = cxEnabler.getMolParam(fus3SpecTr);
	const bnd::modMolState& rFus3State = pFus3->externState(fus3Prm);

	// Generate new reaction only if Fus3 is bound to ATP and is
	// sufficiently phosphorylated
	if((rFus3State[atpModNdx] == pAtpBound)
	   && (minFus3PhosCount <= rFus3State.modCount(pPhosphorylated)))
	  {
	    // Get the state of Ste7, in order to determine which
	    // phosphorylation sites are available to Fus3.
	    plx::plexMolSpec ste7SpecTr = cxEnabler.translateMolSpec(ste7Spec);
	    bnd::molParam ste7Prm = cxEnabler.getMolParam(ste7SpecTr);
	    const bnd::modMolState& rSte7State = pSte7->externState(ste7Prm);

	    // Generate a reaction for each modification site on Ste7
	    // that can be phosphorylated by Fus3 and is currently free.
	    int targetModSiteNdx = activityMask.size();
	    while(0 < targetModSiteNdx--)
	      {
		if(activityMask[targetModSiteNdx]
		   && (rSte7State[targetModSiteNdx] == pNone))
		  {
		    // Construct a phosphorylation reaction.
		    mzr::reaction* pReaction = new mzr::reaction();
		    pFamily->addReaction(pReaction,
					 rMzrUnit);

		    // Add the notifying complex as the only substrate of the
		    // reaction.
		    pReaction->addSubstrate(cxEnabler.getSpecies(),
					    1);
		    pReaction->sensitizeToSubstrates(rMzrUnit);

		    // Set the (unary) reaction rate.
		    pReaction->setRate(pExtrap->getRate(rContext));

		    // Make the new state of Ste7 by phosphoryating the free
		    // modification site.
		    bnd::modMolState nuSte7State(rSte7State);
		    nuSte7State[targetModSiteNdx] = pPhosphorylated;
		    bnd::molParam nuSte7Prm = pSte7->internState(nuSte7State);

		    // Make the new state of Fus3 by replacing bound ATP with
		    // bound ADP.
		    bnd::modMolState nuFus3State(rFus3State);
		    nuFus3State[atpModNdx] = pAdpBound;
		    bnd::molParam nuFus3Prm = pFus3->internState(nuFus3State);

		    // Assemble the molParams of the product complex.
		    std::vector<bnd::molParam>
		      nuMolParams(cxEnabler.getMolParams());
		    nuMolParams[ste7SpecTr] = nuSte7Prm;
		    nuMolParams[fus3SpecTr] = nuFus3Prm;

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

