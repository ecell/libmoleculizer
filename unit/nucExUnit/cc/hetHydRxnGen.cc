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
#include "nucEx/hetHydRxnGen.hh"

namespace nucEx
{
  void
  hetHydRxnGen::makeReactions(const plx::omniInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the new complex in accessors.
	plx::cxOmni cxHetHyd(rContext);

	// Get the index of the nucleotide binding protein mol in the complex
	// of whose appearance we're being notified.  We also need this later
	// to construct the new molParams of the post-reaction complex.
	plx::plexMolSpec targetMolSpec = cxHetHyd.translateMolSpec(nbMolSpec);

	// Get the state of the nucleotide-binding protein.  To get the
	// nucleotide binding mol, we must first translate the target mol's
	// index using the injection of the "enabling subcomplex" of the
	// hydrolysing mol and the "assistant" mol.
	bnd::molParam hydMolParam = cxHetHyd.getMolParam(targetMolSpec);
	const bnd::modMolState& rHydMolState = pMol->externState(hydMolParam);

	// Generate a reaction only if the mol is bound to the unhydrolysed
	// nucleotide.
	if(rHydMolState[modSiteNdx] == pUnhydrolysedBound)
	  {
	    // Create the reaction.
	    mzr::reaction* pReaction = new mzr::reaction();
	    pFamily->addReaction(pReaction,
				 rMzrUnit);

	    // Add the complex containing the nucleotide-binding protein and
	    // the "assitant" protein as the only substrate.
	    pReaction->addSubstrate(cxHetHyd.getSpecies(),
				    1);
	    pReaction->sensitizeToSubstrates(rMzrUnit);

	    // Set the reaction rate.  (This is basically a unary reaction.)
	    pReaction->setRate(pExtrap->getRate(rContext));

	    // Add phosphate as a product of the reaction.
	    pReaction->addProduct(pPhosphate,
				  1);

	    // Make the new state of the nucleotide-binding protein, now bound
	    // to the hydrolysed nucleotide.
	    bnd::modMolState nuState(rHydMolState);
	    nuState[modSiteNdx] = pHydrolysedBound;
	    bnd::molParam nuParam = pMol->internState(nuState);

	    // Make parameter for the post-reaction complex.
	    std::vector<bnd::molParam> molParams(cxHetHyd.getMolParams());
	    molParams[targetMolSpec] = nuParam;

	    // Make the post-reaction complex.
	    plx::plexSpecies* pPost
	      = cxHetHyd.getPlexFamily().getMember(molParams);

	    // Add it as a product of multiplicity 1.
	    pReaction->addProduct(pPost,
				  1);
	  }
      }
  }
}
