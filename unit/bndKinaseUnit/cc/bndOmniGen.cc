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
#include "mol/modMol.hh"
#include "plex/plexUnit.hh"
#include "bndKinase/bndOmniGen.hh"
#include "bndKinase/bndKinaseXcpt.hh"

namespace bndKinase
{
  class exchangeSmallMols :
    public std::unary_function<smallMolExchange, void>
  {
    const plx::cxOmni& rEnabling;
    std::vector<bnd::mol*>& rMols;
  public:
    exchangeSmallMols(const plx::cxOmni& rEnablingOmni,
		      std::vector<bnd::mol*>& rResultMols) :
      rEnabling(rEnablingOmni),
      rMols(rResultMols)
    {}

    void
    operator()(const smallMolExchange& rExchange) const
    {
      // Get the index, in the actual reactant species, of the mol to be
      // exchanged.
      plx::plexMolSpec exchangedMolSpecTr
	= rEnabling.translateMolSpec(rExchange.exchangedMolSpec);

      // Swap in the new small-mol.
      rMols[exchangedMolSpecTr] = rExchange.pReplacementMol;
    }
  };

  class exchangeSmallMolParams :
    public std::unary_function<smallMolExchange, void>
  {
    const plx::cxOmni& rEnabling;
    const plx::plexIsoPair& rResultToResultParadigm;
    std::vector<bnd::molParam>& rProductMolParams;
  public:
    exchangeSmallMolParams(const plx::cxOmni& rEnablingOmni,
			   const plx::plexIsoPair& refResultToResultParadigm,
			   std::vector<bnd::molParam>& refProductMolParams) :
      rEnabling(rEnablingOmni),
      rResultToResultParadigm(refResultToResultParadigm),
      rProductMolParams(refProductMolParams)
    {}

    void
    operator()(const smallMolExchange& rExchange) const
    {
      // Translate the molSpec of the exchanged mol into "result" indexing.
      plx::plexMolSpec exchangedMolSpecTr
	= rEnabling.translateMolSpec(rExchange.exchangedMolSpec);

      // Determine where the default state of the swapped-in
      // small-mol should go in "paradigm" indexing.
      plx::plexMolSpec paradigmSmallMolSpec
	= rResultToResultParadigm.forward.molMap[exchangedMolSpecTr];

      // Install the default state of the swapped in small-mol.
      rProductMolParams[paradigmSmallMolSpec]
	= rExchange.pReplacementMol->getDefaultParam();
    }
  };

  class exchangeModifications :
    public std::unary_function<modificationExchange, void>
  {
    const plx::plexIsoPair& rResultToResultParadigm;
    const std::vector<bnd::mol*>& rResultParadigmMols;
    std::vector<bnd::molParam>& rProductMolParams;
  public:
    exchangeModifications(const plx::plexIsoPair& refResultToResultParadigm,
			  const std::vector<bnd::mol*>& refResultParadigmMols,
			  std::vector<bnd::molParam>& refProductMolParams) :
      rResultToResultParadigm(refResultToResultParadigm),
      rResultParadigmMols(refResultParadigmMols),
      rProductMolParams(refProductMolParams)
    {}

    void
    operator()(const modificationExchange& rExchange) const
    {
      // Determine what mol the modification will be swapped in.
      plx::plexMolSpec paradigmModMolSpec
	= rResultToResultParadigm.forward.molMap[rExchange.modMolSpec];

      // Get the mol and make sure it's a mod-mol.
      bnd::modMol* pModMol
	= dynamic_cast<bnd::modMol*>(rResultParadigmMols[paradigmModMolSpec]);

      if(! pModMol)
	throw exchangedModXcpt(rResultParadigmMols[paradigmModMolSpec]);

      // Get the current modification state.
      bnd::modMolState exchangedMolState
	= pModMol->externState(rProductMolParams[paradigmModMolSpec]);

      // Substitute in the modification.
      exchangedMolState[rExchange.modSiteNdx] = rExchange.pReplacementMod;

      // Intern the new modMol state.
      const bnd::molParam nuMolParam
	= pModMol->internState(exchangedMolState);

      // Install the molParam into the vector of molParams for the primary
      // product species.
      rProductMolParams[paradigmModMolSpec] = nuMolParam;
    }
  };
  
  void
  bndOmniRxnGen::
  makeReactions(const plx::omniInContext& rContext) const
  {
    // See if reaction generation is enabled.
    if(rMzrUnit.generateOk)
      {
	// Wrap the enabling complex in accessors.
	plx::cxOmni cxEnabling(rContext);

	// Construct the reaction and intern it for memory management.
	mzr::reaction* pReaction = new mzr::reaction();
	pFamily->addReaction(pReaction,
			     rMzrUnit);

	// Add the triggering complex as a reactant of multiplicity 1.
	pReaction->addSubstrate(cxEnabling.getSpecies(),
				1);

	// If an auxiliary reactant was specified, add it with multiplicity 1.
	if(pAdditionalReactant)
	  {
	    pReaction->addSubstrate(pAdditionalReactant,
				    1);
	  }

	// Sensitize the reaction to its reactants (i.e. add the reaction to
	// the sensitivity list of each reaction, and to the sensitivity list
	// of a few other state variables, such as volume.)
	pReaction->sensitizeToSubstrates(rMzrUnit);

	// If an auxiliary product was specified, add it with multiplicity 1.
	if(pAdditionalProduct)
	  {
	    pReaction->addProduct(pAdditionalProduct,
				  1);
	  }

	// Set the rate of the reaction.
	pReaction->setRate(pExtrapolator->getRate(rContext));

	// Now construct the primary product species, the modified complex.
	// This occurs in two steps, making the structural complex and making
	// the parameters for the mols in the complex.

	// First, we have to construct its vector of mols, which is obtained
	// from the primary reactant's vector of mols by performing the
	// small-mol exchanges, if any.
	plx::plex resultPlex(cxEnabling.getPlexFamily().getParadigm());

	std::for_each(smallMolExchanges.begin(),
		      smallMolExchanges.end(),
		      exchangeSmallMols(cxEnabling,
					resultPlex.mols));

	// Recognize the product complex, retaining the recognition map so
	// that we can reorder the states of the mols appropriately.
	plx::plexIsoPair resultToResultParadigm;
	plx::plexFamily* pProductFamily
	  = rPlexUnit.recognize(resultPlex,
				resultToResultParadigm);

	// Now make the parameters for the mols in the product complex.

	// We start with a copy of the triggering complex's mol parameters.
	const std::vector<bnd::molParam>& rEnablingMolParams
	  = cxEnabling.getMolParams();
	int paradigmMolNdx = rEnablingMolParams.size();
	std::vector<bnd::molParam> resultMolParams(paradigmMolNdx);
	while(0 < paradigmMolNdx--)
	  {
	    // Get the index of the mol in the triggering reactant complex.
	    int molNdx
	      = resultToResultParadigm.backward.molMap[paradigmMolNdx];

	    resultMolParams[paradigmMolNdx] = rEnablingMolParams[molNdx];
	  }

	// Substitute the default parameters for all small-mols that were
	// substituted in.
	std::for_each(smallMolExchanges.begin(),
		      smallMolExchanges.end(),
		      exchangeSmallMolParams(cxEnabling,
					     resultToResultParadigm,
					     resultMolParams));

	// Perform the modification exchanges.
	std::for_each(modificationExchanges.begin(),
		      modificationExchanges.end(),
		      exchangeModifications(resultToResultParadigm,
					    pProductFamily->getParadigm().mols,
					    resultMolParams));

	// Construct the primary product species.
	plx::plexSpecies* pProductSpecies
	  = pProductFamily->getMember(resultMolParams);

	// Add the primary product species to the reaction, with multiplicity
	// 1.
	pReaction->addProduct(pProductSpecies,
			      1);

      }
  }
}
