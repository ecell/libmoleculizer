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
#include "modKinase/nucBindRxnGen.hh"

namespace kinase
{
  void
  nucBindRxnGen::makeReactions(const plx::molInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the new context in accessors.
	plx::cxMol cxBinder(rContext);

	// Determine if the mol is bound to nucleotide, in which case we
	// generate an unbinding reaction, bound to nothing, in which case we
	// generate a binding reaction, or bound to xyz, in which case we do
	// nothing.
	//
	// Here, bound to xxx means that the appropriate modification is at the
	// modification site that is being used to model nucleotide binding.
	bnd::molParam binderParam = cxBinder.getMolParam();
	const bnd::modMolState& rBinderState
	  = pMol->externState(binderParam);

	const bnd::modification* pCurrentMod = rBinderState[modSiteNdx];
	if(pCurrentMod == pBound)
	  {
	    // Generate an unbinding reaction.
	    mzr::reaction* pUnbinding = new mzr::reaction();
	    pFamily->addReaction(pUnbinding,
				 rMzrUnit);

	    // The complex containing the nucleotide binding mol is the
	    // only substrate.
	    pUnbinding->addSubstrate(cxBinder.getSpecies(),
				     1);
	    pUnbinding->sensitizeToSubstrates(rMzrUnit);

	    // Extrapolate the (unary) reaction rate.
	    pUnbinding->setRate(pNucUnbindExtrap->getRate(rContext));

	    // Add nucleotide as a product.
	    pUnbinding->addProduct(pNuc,
				   1);

	    // Create the unbound state of the nucleotide-binding mol.
	    bnd::modMolState nuState(rBinderState);
	    nuState[modSiteNdx] = pNone;
	    bnd::molParam nuParam = pMol->internState(nuState);

	    // Create the new state of the complex.
	    std::vector<bnd::molParam> molParams(cxBinder.getMolParams());
	    molParams[cxBinder.getMolSpec()] = nuParam;

	    // Make the product species.
	    plx::plexSpecies* pResult
	      = cxBinder.getPlexFamily().getMember(molParams);

	    // Add it as a product of multiplicity 1.
	    pUnbinding->addProduct(pResult,
				   1);
	  }
	else if(pCurrentMod == pNone)
	  {
	    // Generate a binding reaction.
	    mzr::reaction* pBinding = new mzr::reaction();
	    pFamily->addReaction(pBinding,
				 rMzrUnit);

	    // Add the substrates.
	    pBinding->addSubstrate(cxBinder.getSpecies(),
				   1);
	    pBinding->addSubstrate(pNuc,
				   1);
	    pBinding->sensitizeToSubstrates(rMzrUnit);

	    // Extrapolate the reaction rate.
	    pBinding->setRate(pNucBindExtrap->getRate(rContext));

	    // Create the bound state of the nucleotide-binding mol.
	    bnd::modMolState nuState(rBinderState);
	    nuState[modSiteNdx] = pBound;
	    bnd::molParam nuParam = pMol->internState(nuState);

	    // Create the new state of the complex.
	    std::vector<bnd::molParam> molParams(cxBinder.getMolParams());
	    molParams[cxBinder.getMolSpec()] = nuParam;

	    // Make the product species.
	    plx::plexSpecies* pResult
	      = cxBinder.getPlexFamily().getMember(molParams);

	    // Add it as a product of multiplicity 1.
	    pBinding->addProduct(pResult,
				 1);
	  }
      }
  }
}
