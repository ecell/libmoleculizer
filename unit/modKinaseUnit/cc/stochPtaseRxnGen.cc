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
#include "modKinase/stochPtaseRxnGen.hh"

namespace kinase
{
  void
  stochPtaseRxnGen::makeReactions(const plx::molInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	
	// Wrap the new context in accessors.
	plx::cxMol cxSubstrate(rContext);

	// Get the state of the mol.
	bnd::molParam substrateParam = cxSubstrate.getMolParam();
	const bnd::modMolState& rSubstrateState = pMol->externState(substrateParam);

	// Run through the activity mask and modification sites, generating
	// a reaction for each phosphorylated modification site.
	//
	// Once again, I need a binary version of for_each.
	for(int modSiteNdx = 0;
	    modSiteNdx < pMol->modSiteCount();
	    ++modSiteNdx)
	  {
	    if(phosMask[modSiteNdx]
	       && (rSubstrateState[modSiteNdx] == pPhosMod))
	      {
		// Create the new reaction and install it in the reaction family
		// for memory management.
		mzr::reaction* pReaction = new mzr::reaction();
		pFamily->addReaction(pReaction,
				     rMzrUnit);

		// Add the substrate-containing complex and the phosphatase
		// as substrates of the reaction.
		pReaction->addSubstrate(cxSubstrate.getSpecies(),
					1);
		pReaction->addSubstrate(pPtase,
					1);
		pReaction->sensitizeToSubstrates(rMzrUnit);

		// Extrapolate the (binary) reaction rate.
		pReaction->setRate(pExtrap->getRate(rContext));

		// Add phosphate as a product.
		pReaction->addProduct(pPhosphate,
				      1);

		// Make the post-reaction state of the substrate protein, with
		// the phosphate removed from the modSiteNdx-th modification
		// site.
		bnd::modMolState nuState(rSubstrateState);
		nuState[modSiteNdx] = pNoneMod;
		bnd::molParam nuParam = pMol->internState(nuState);

		// Make parameters for the post-reaction complex.
		std::vector<bnd::molParam>
		  molParams(cxSubstrate.getMolParams());
		molParams[cxSubstrate.getMolSpec()] = nuParam;

		// Make the post reaction complex species.
		plx::plexSpecies* pPost
		  = cxSubstrate.getPlexFamily().getMember(molParams);

		// Add it as a product of multiplicity 1.
		pReaction->addProduct(pPost,
				      1);
	      }
	  }
      }
  }
}
