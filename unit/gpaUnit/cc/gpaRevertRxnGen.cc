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
#include "gpa/gpaRevertRxnGen.hh"

namespace gpa
{
  void
  gpaRevertRxnGen::makeReactions(const plx::molInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the new context in accessors.
	plx::cxMol cxReverter(rContext);

	// Get the state of the gpa mol.
	bnd::molParam pState = cxReverter.getMolParam();
	const bnd::modMolState& rModMolState = pMol->externState(pState);
  
	// Generate a reaction only if the mol is GTP-bound.
	if(rModMolState[gtpSiteNdx] == pGtpBound)
	  {
	    // Create the new reaction and install it in the reaction family
	    // for memory management.
	    mzr::reaction* pReaction = new mzr::reaction();
	    pFamily->addReaction(pReaction,
				 rMzrUnit);

	    // Add the complex containing the gpa protein as a substrate.
	    pReaction->addSubstrate(cxReverter.getSpecies(),
				    1);
	    pReaction->sensitizeToSubstrates(rMzrUnit);

	    // Set the reaction rate.
	    pReaction->setRate(pExtrap->getRate(rContext));

	    // Add phosphate as a product.
	    pReaction->addProduct(pPhosphate,
				  1);

	    // Make the new (gdp) state of the gpa mol.
	    bnd::modMolState nuState(rModMolState);
	    nuState[gtpSiteNdx] = pGdpBound;
	    bnd::molParam nuParam = pMol->internState(nuState);

	    // Make the product species param.
	    std::vector<bnd::molParam> molParams(cxReverter.getMolParams());
	    molParams[cxReverter.getMolSpec()] = nuParam;

	    // Make the product species.
	    plx::plexSpecies* pResult
	      = cxReverter.getPlexFamily().getMember(molParams);

	    // Add it as a product of multiplicity 1.
	    pReaction->addProduct(pResult,
				  1);
	  }
      }
  }
}
