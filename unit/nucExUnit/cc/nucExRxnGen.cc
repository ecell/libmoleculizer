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
#include "nucEx/nucExRxnGen.hh"

namespace nucEx
{
  // Tests to see if a plexIsoPair carrys one given plexMolSpec to another.
  class imbeddingHitsTargetMol :
    public std::unary_function<plx::subPlexSpec, bool>
  {
    plx::plexMolSpec enablingSpec;
    plx::plexMolSpec actualSpec;
  public:
    imbeddingHitsTargetMol(plx::plexMolSpec enablingMolSpec,
			   plx::plexMolSpec actualMolSpec) :
      enablingSpec(enablingMolSpec),
      actualSpec(actualMolSpec)
    {}

    bool
    operator()(const argument_type& rSubPlexSpec) const
    {
      const plx::plexIsoPair& rInjection = rSubPlexSpec.injection();
      return rInjection.forward.applyToMolSpec(enablingSpec) == actualSpec;
    }
  };

  void
  nucExRxnGen::makeReactions(const plx::molInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	
	plx::cxMol cxExchangeMol(rContext);

	// Determine whether we are in the "enabled" situation or not.
	// 
	// See if the target mol is situated appropriately in the "enabling"
	// complex.  If so, then one set of reaction rates is used, if not,
	// then another.  The complex might more appropriately be called the
	// "enhancing" complex.
	plx::plexFamily& rNotifyingFamily = cxExchangeMol.getPlexFamily();
	plx::plexMolSpec exchangeMolSpec = cxExchangeMol.getMolSpec();

	// This should probably change under the new omniPlex regime.
	// I've fixed "getOmniImbeddings" so that it compiles, but not
	// worked on the solution to this...
	zzzzzz;
	
	std::vector<plx::subPlexSpec> imbeddings
	  = rNotifyingFamily.getOmniImbeddings(pEnablingFamily);

	std::vector<plx::subPlexSpec>::iterator iAppropriate
	  = std::find_if(imbeddings.begin(),
			 imbeddings.end(),
			 imbeddingHitsTargetMol(enablingSpec,
						exchangeMolSpec));

	bool enabled = (imbeddings.end() != iAppropriate);

	// Look at the state of the target mol.  Is it right to generate
	// a binding reaction, an unbinding reaction, or neither?
	bnd::molParam exchangeMolParam = cxExchangeMol.getMolParam();
	const bnd::modMolState& rExchangeMolState
	  = pMol->externState(exchangeMolParam);

	// If the target mol is unmodified at the "nucleotide-binding"
	// modification site, then we can generate a binding reaction.
	const bnd::modification* pCurrentMod = rExchangeMolState[modSiteNdx];

	if(pCurrentMod == pNoneMod)
	  {
	    // Create nucleotide binding reaction.
	    mzr::reaction *pReaction = new mzr::reaction();
	    pFamily->addReaction(pReaction,
				 rMzrUnit);

	    // Add its substrates.
	    pReaction->addSubstrate(cxExchangeMol.getSpecies(),
				    1);
	    pReaction->addSubstrate(pNucSpecies,
				    1);
	    pReaction->sensitizeToSubstrates(rMzrUnit);

	    // Set its rate, which depends on whether we're "enabled."
	    // 	double onConst = enabled ? enabledOnConst : plainOnConst;
	    // 	pReaction->setRate(mzr::bindingRate(onConst,
	    // 				       cxExchangeMol.getPlexWeight(),
	    // 				       pNucSpecies->getWeight()));

	    // Extrapolate the rate of binding, which depends on whether
	    // we're in the enabling complex or not.
	    nucExBindExtrapolator* pExtrap
	      = enabled ? pEnabledBindExtrap : pPlainBindExtrap;

	    pReaction->setRate(pExtrap->getRate(rContext));

	    // Make the new state of the nucleotide binding mol.
	    bnd::modMolState nuState(rExchangeMolState);
	    nuState[modSiteNdx] = pBoundMod;
	    bnd::molParam nuParam = pMol->internState(nuState);

	    // Make the product species parameters.
	    std::vector<bnd::molParam> molParams(cxExchangeMol.getMolParams());
	    molParams[exchangeMolSpec] = nuParam;

	    // Make the product species.
	    plx::plexSpecies* pResult
	      = rNotifyingFamily.getMember(molParams);

	    // Add it as a product of multiplicity 1.
	    pReaction->addProduct(pResult,
				  1);
	  }

	// If the target mol is already bound to nucleotide, construct an
	// unbinding reaction.
	else if(pCurrentMod == pBoundMod)
	  {
	    // Create nucleotide unbinding reaction.
	    mzr::reaction* pReaction = new mzr::reaction();
	    pFamily->addReaction(pReaction,
				 rMzrUnit);

	    // Add its substrate.
	    pReaction->addSubstrate(cxExchangeMol.getSpecies(),
				    1);
	    pReaction->sensitizeToSubstrates(rMzrUnit);

	    // Set its rate, which depends on whether we're "enabled."
	    // 	double offRate = enabled ? enabledOffRate : plainOffRate;
	    // 	pReaction->setRate(offRate);

	    // Extrapolate the rate of unbinding, which depends on whether
	    // we're in the enabling complex or not.
	    nucExUnbindExtrapolator* pExtrap
	      = enabled ? pEnabledUnbindExtrap : pPlainUnbindExtrap;
	
	    pReaction->setRate(pExtrap->getRate(rContext));

	    // Add the nucleotide as a product.
	    pReaction->addProduct(pNucSpecies,
				  1);

	    // Make the new unbound state of the nucleotide binding mol.
	    bnd::modMolState nuState(rExchangeMolState);
	    nuState[modSiteNdx] = pNoneMod;
	    bnd::molParam nuParam = pMol->internState(nuState);

	    // Make product species parameters for nucleotide-unbound version
	    // of notifying complex.
	    std::vector<bnd::molParam> molParams(cxExchangeMol.getMolParams());
	    molParams[exchangeMolSpec] = nuParam;

	    // Make the product species.
	    plx::plexSpecies* pResult
	      = rNotifyingFamily.getMember(molParams);

	    // Add it as a product of multiplicity 1.
	    pReaction->addProduct(pResult,
				  1);
	  }

	// Otherwise, the "nucleotide-binding" modification site has some
	// unexpected modification.  Issue a warning?
      }
  }
}
