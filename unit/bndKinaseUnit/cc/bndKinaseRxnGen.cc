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
#include "bndKinase/bndKinaseRxnGen.hh"

namespace bndKinase
{
  void
  bndKinaseRxnGen::makeReactions(const plx::omniInContext& rContext) const
  {
    // See if reaction generation is turned off.
    if(rMzrUnit.generateOk)
      {
	// Wrap the enabling complex in accessors.
	plx::cxOmni cxEnabling(rContext);

	// Construct the reation and intern it for memory management.
	mzr::reaction* pReaction = new mzr::reaction();
	pFamily->addReaction(pReaction,
			     rMzrUnit);

	// Add the complex as the only substrate of the reaction, with
	// multiplicity 1.
	pReaction->addSubstrate(cxEnabling.getSpecies(),
				1);

	// (Note that mzrUnit is used here to sensitize to volume, etc.)
	pReaction->sensitizeToSubstrates(rMzrUnit);

	// Set the (unary) reaction rate.
	pReaction->setRate(pExtrap->getRate(rContext));

	// Construct the product complex.  Its mols are the same as the
	// mols of the triggering complex, except that ATP is replaced
	// with ADP.  Its bindings are the same as the bindings of the
	// triggering complex.
	plx::plex resultPlex(cxEnabling.getPlexFamily().getParadigm());
	plx::plexMolSpec atpMolSpecTr = cxEnabling.translateMolSpec(atpSpec);
	resultPlex.mols[atpMolSpecTr] = pADPmol;

	// Recognize the product complex, retaining the recognition map
	// so that we can reorder the states of the mols appropriately......
	plx::plexIsoPair resultToParadigm;
	plx::plexFamily* pProductFamily
	  = rPlexUnit.recognize(resultPlex,
				resultToParadigm);

	// Construct the vector of molParams by permuting the mols of the
	// reactant species with the recogntion map, then substituting in
	// the new state of the substrate mol.
	int paradigmMolNdx = cxEnabling.getMolParams().size();
	std::vector<bnd::molParam> productMolParams(paradigmMolNdx);
	while(0 < paradigmMolNdx--)
	  {
	    // Get the index of the mol in the triggering complex
	    // that corresponds to the mol at paradigmMolNdx
	    // in the paradigm of the product complex.
	    int molNdx = resultToParadigm.backward.molMap[paradigmMolNdx];

	    // Copy the molparams from the reactant species into their
	    // corresponding posistions for the result species.
	    productMolParams[paradigmMolNdx]
	      = cxEnabling.getMolParams()[molNdx];
	  }

	// Get the state of the substrate mol.
	plx::plexMolSpec substrateSpecTr
	  = cxEnabling.translateMolSpec(substrateSpec);
	bnd::molParam substrateMolParam
	  = cxEnabling.getMolParam(substrateSpecTr);
	const bnd::modMolState& rSubstrateState
	  = pSubstrate->externState(substrateMolParam);

	// Construct the new, phosphorylated state of the substrate mol.
	bnd::modMolState nuSubstrateState(rSubstrateState);
	nuSubstrateState[phosSiteNdx] = pPhosphorylated;
	bnd::molParam nuSubstrateParam
	  = pSubstrate->internState(nuSubstrateState);

	// Substitute the new state of the substrate mol into the
	// productMolParams.
	//
	// The index of the substrate mol in the constructed result
	// complex is substrateSpecTr, the same as in the reactant
	// species.  We need to put the new parameter for this mol in the
	// paradigm of the result species at the index associated to
	// substrateSpecTr by the recognition map resultToParadigm.
	plx::plexMolSpec paradigmSubstrateSpec
	  = resultToParadigm.forward.molMap[substrateSpecTr];
	productMolParams[paradigmSubstrateSpec] = nuSubstrateParam;

	// Substitute the default state of the ADP small-mol for the
	// state of the ATP small-mol.
	plx::plexMolSpec paradigmADPSpec
	  = resultToParadigm.forward.molMap[atpMolSpecTr];
	productMolParams[paradigmADPSpec] = pADPmol->getDefaultState();

	// Construct the product species.
	plx::plexSpecies* pProductSpecies
	  = pProductFamily->getMember(productMolParams);

	// Add the product species to the reaction with multiplicity 1.
	pReaction->addProduct(pProductSpecies,
			      1);
      }
  }
}
