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

#include "mol/modMol.hh"
#include "plex/cxOmniParam.hh"
#include "plex/plexUnit.hh"
#include "bndKinase/modRxnGen.hh"

namespace bndKinase
{
  void
  modRxnGen::makeReactions(const plx::omniInContext& rContext) const
  {
    // See if reaction generation is turned off.
    if(rMzrUnit.generateOk)
      {
	// Wrap the enabling omni in accessors.
	plx::cxOmni cxEnabling(rContext);

	// Construct the generic modification reaction and intern it
	// for memory management.
	mzr::reaction* pReaction = new mzr::reaction();
	pFamily->addReaction(pReaction,
			     rMzrUnit);

	// Add the complex as the only reactant, with multiplicity 1.
	pReaction->addSubstrate(cxEnabling.getSpecies(),
				1);

	// Sensitize to reactants, volume, etc.
	pReaction->sensitizeToSubstrates(rMzrUnit);

	// Set the extrapolated (unary) reaction rate.
	pReaction->setRate(pExtrap->getRate(rContext));

	// Make a copy of the molParams of the notifying species.
	std::vector<bnd::molParam>
	  productMolParams(cxEnabling.getMolParams());

	// Get the state of the substrate mol.
	plx::plexMolSpec substrateSpecTr
	  = cxEnabling.translateMolSpec(substrateSpec);
	bnd::molParam substrateMolParam
	  = cxEnabling.getMolParam(substrateSpecTr);
	const bnd::modMolState& rSubstrateState
	  = pSubstrate->externState(substrateMolParam);

	// Construct the new modification state of the susbtrate mol
	// by resetting the modification state at the designated site.
	bnd::modMolState nuSubstrateState(rSubstrateState);
	nuSubstrateState[substrateModSiteNdx] = pInstalledMod;
	bnd::molParam nuSubstrateParam
	  = pSubstrate->internState(nuSubstrateState);

	// Substitute the new state of the substrate mod-mol into the
	// productMolParams.
	productMolParams[substrateSpecTr] = nuSubstrateParam;

	// Construct the product species.
	plx::plexSpecies* pProductSpecies
	  = cxEnabling.getPlexFamily().getMember(productMolParams);

	// Add the product species to the reaction with multiplicity 1.
	pReaction->addProduct(pProductSpecies,
			      1);
      }
  }
}
