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

#include "mzr/pchem.hh"
#include "mzr/mzrUnit.hh"
#include "plex/cxMolParam.hh"
#include "modKinase/modKinaseRxnGen.hh"
#include "modKinase/kinaseXcpt.hh"

namespace kinase
{
  void modKinaseRxnGenPair::
  makeBinaryReactions(const plx::molInContext& rKinaseContext,
		      const plx::molInContext& rSubstrateContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the mol contexts in accessor functions.
	plx::cxMol cxKinase(rKinaseContext);
	plx::cxMol cxSubstrate(rSubstrateContext);

	// Get the modification state of the kinase.  Later, if all mols
	// become modMols, this dynamic_cast will go away.
	bnd::molParam kinaseParam = cxKinase.getMolParam();
	const bnd::modMolState& rKinaseState
	  = pKinase->externState(kinaseParam);

	// We only generate a reaction if the kinase is bound to ATP.
	if(rKinaseState[atpSiteNdx] == pATPbound)
	  {
	    // Get the modification state of the substrate.  Later, if
	    // all mols become modMols, this dynamic_cast can go away.
	    bnd::molParam substrateParam = cxSubstrate.getMolParam();
	    const bnd::modMolState& rSubstrateState
	      = pSubstrate->externState(substrateParam);

	    // For debugging: check that the number of modification
	    // sites on the substrate is the same as the length of
	    // the activity mask.
	    if(rSubstrateState.size() != activityMask.size())
	      throw maskModSiteCountXcpt();

	    // Run down the mask.
	    for(int modSiteNdx = 0;
		modSiteNdx < (int) activityMask.size();
		++modSiteNdx)
	      {
		// Is the (phosphorylation) site "enabled" and is it
		// free?
		if(activityMask[modSiteNdx]
		   && (rSubstrateState[modSiteNdx] == pNone))
		  {
		    // Create a new reaction and add it to the family.
		    mzr::reaction* pReaction = new mzr::reaction();
		    pFamily->addReaction(pReaction,
					 rMzrUnit);

		    // Convert the weight-independent reactionConst into the
		    // reaction rate appropriate for (the weights of) these
		    // reactants.  double rate =
		    // mzr::bindingRate(reactionConst,
		    // cxKinase.getPlexWeight(), cxSubstrate.getPlexWeight());

		    pReaction->setRate(pExtrap->getRate(rKinaseContext,
							rSubstrateContext));

		    // Add the kinase and the substrate as substrates(!)
		    // and sensitize the new reaction to them.
		    pReaction->addSubstrate(cxKinase.getSpecies(),
					    1);
		    pReaction->addSubstrate(cxSubstrate.getSpecies(),
					    1);
		    pReaction->sensitizeToSubstrates(rMzrUnit);

		    // Construct the ADP form of the kinase's parameter.
		    // This should be done outside this loop.
		    bnd::modMolState nuKinaseState(rKinaseState);
		    nuKinaseState[atpSiteNdx] = pADPbound;
		    bnd::molParam nuKinaseParam
		      = pKinase->internState(nuKinaseState);

		    // Make the product species generator for the ADP form
		    // of the kinase complex, and install it in the reaction.
		    std::vector<bnd::molParam>
		      kinaseOutParams(cxKinase.getMolParams());
		    kinaseOutParams[cxKinase.getMolSpec()] = nuKinaseParam;

		    // Construct the ADP form of the kinase complex.
		    plx::plexSpecies* pADPKinase
		      = cxKinase.getPlexFamily().getMember(kinaseOutParams);

		    // Add it as a product of multiplicity 1.
		    pReaction->addProduct(pADPKinase,
					  1);

		    // Construct the phosphorylated form of the substrate
		    // mol param.
		    bnd::modMolState nuSubstrateState(rSubstrateState);
		    nuSubstrateState[modSiteNdx] = pPhosphorylated;
		    bnd::molParam nuSubstrateParam
		      = pSubstrate->internState(nuSubstrateState);

		    // Make parameters of phosphorylated substrate complex.
		    std::vector<bnd::molParam>
		      substrateOutParams(cxSubstrate.getMolParams());
		    substrateOutParams[cxSubstrate.getMolSpec()]
		      = nuSubstrateParam;

		    // Create the phosphorylated substrate complex.
		    plx::plexSpecies* pSubstrateOut
		      = cxSubstrate.getPlexFamily().getMember(substrateOutParams);

		    // Add it as a product of multiplicity 1.
		    pReaction->addProduct(pSubstrateOut,
					  1);
		  }
	      }
	  }
      }
  }
}
