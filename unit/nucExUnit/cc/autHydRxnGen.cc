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
#include "plex/cxMolParam.hh"
#include "nucEx/autHydRxnGen.hh"

namespace nucEx
{
  void
  autHydRxnGen::makeReactions(const plx::molInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the new context in accessors.
	plx::cxMol cxHydrolyser(rContext);

	// Determine if the mol is bound to the unhydrolysed nucleotide
	// (GTP).
	bnd::molParam hydrolyserParam = cxHydrolyser.getMolParam();
	const bnd::modMolState& rHydrolyserState
	  = pMol->externState(hydrolyserParam);

	if(rHydrolyserState[modSiteNdx] == pUnhydrolysedBound)
	  {
	    // Create the new reaction ans install it in it reaction family
	    // for memory management.
	    mzr::reaction* pReaction = new mzr::reaction();
	    pFamily->addReaction(pReaction,
				 rMzrUnit);

	    // Add the complex containing the nuclotide binding protein
	    // as the only substrate.
	    pReaction->addSubstrate(cxHydrolyser.getSpecies(),
				    1);
	    pReaction->sensitizeToSubstrates(rMzrUnit);

	    // Set the reaction rate.
	    pReaction->setRate(pExtrap->getRate(rContext));

	    // Add phosphate as a product.
	    pReaction->addProduct(pPhosphate,
				  1);

	    // Make the post-reaction state of the nucleotide-binding protein.
	    bnd::modMolState nuState(rHydrolyserState);
	    nuState[modSiteNdx] = pHydrolysedBound;
	    bnd::molParam nuParam = pMol->internState(nuState);

	    // Make parameter for the post-reaction complex.
	    std::vector<bnd::molParam> molParams(cxHydrolyser.getMolParams());
	    molParams[cxHydrolyser.getMolSpec()] = nuParam;

	    // Make the post reaction complex.
	    plx::plexSpecies* pPost
	      = cxHydrolyser.getPlexFamily().getMember(molParams);

	    // Add it as a product of multiplicity 1.
	    pReaction->addProduct(pPost,
				  1);
	  }
      }
  }
}
