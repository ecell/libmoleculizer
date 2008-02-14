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

#include "cpt/cptUnit.hh"
#include "cml/cptModMol.hh"
#include "clx/clxUnit.hh"
#include "cft/omniGen.hh"
#include "cft/exchangedModXcpt.hh"

namespace cft
{
  class exchangeSmallMols :
    public std::unary_function<smallMolExchange, void>
  {
    const cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex>& rEnabling;

    std::vector<cml::cptMol*>& rMols;

  public:
    exchangeSmallMols(const cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex>& rEnablingOmni,
		      std::vector<cml::cptMol*>& rResultMols) :
      rEnabling(rEnablingOmni),
      rMols(rResultMols)
    {}

    void
    operator()(const smallMolExchange& rExchange) const
    {
      // Get the index, in the actual reactant species, of the mol to be
      // exchanged.
      cpx::molSpec exchangedMolSpecTr
	= rEnabling.translateMolSpec(rExchange.exchangedMolSpec);

      // Swap in the new small-mol.
      rMols[exchangedMolSpecTr] = rExchange.pReplacementMol;
    }
  };

  class exchangeSmallMolParams :
    public std::unary_function<smallMolExchange, void>
  {
    const cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex>& rEnabling;
    const cpx::plexIso& rResultToResultParadigm;
    std::vector<cpx::molParam>& rProductMolParams;
  public:
    exchangeSmallMolParams(const cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex>& rEnablingOmni,
			   const cpx::plexIso& refResultToResultParadigm,
			   std::vector<cpx::molParam>& refProductMolParams) :
      rEnabling(rEnablingOmni),
      rResultToResultParadigm(refResultToResultParadigm),
      rProductMolParams(refProductMolParams)
    {}

    void
    operator()(const smallMolExchange& rExchange) const
    {
      // Translate the molSpec of the exchanged mol into "result" indexing.
      cpx::molSpec exchangedMolSpecTr
	= rEnabling.translateMolSpec(rExchange.exchangedMolSpec);

      // Determine where the default state of the swapped-in
      // small-mol should go in "paradigm" indexing.
      cpx::molSpec paradigmSmallMolSpec
	= rResultToResultParadigm.forward.molMap[exchangedMolSpecTr];

      // Install the default state of the swapped in small-mol.
      rProductMolParams[paradigmSmallMolSpec]
	= rExchange.pReplacementMol->getDefaultParam();
    }
  };

  class exchangeModifications :
    public std::unary_function<modificationExchange, void>
  {
    const cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex>& rEnabling;
    const cpx::plexIso& rResultToResultParadigm;
    const std::vector<cml::cptMol*>& rResultParadigmMols;
    std::vector<cpx::molParam>& rProductMolParams;
  public:
    exchangeModifications(const cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex>& rEnablingOmni,
			  const cpx::plexIso& refResultToResultParadigm,
			  const std::vector<cml::cptMol*>& refResultParadigmMols,
			  std::vector<cpx::molParam>& refProductMolParams) :
      rEnabling(rEnablingOmni),
      rResultToResultParadigm(refResultToResultParadigm),
      rResultParadigmMols(refResultParadigmMols),
      rProductMolParams(refProductMolParams)
    {}

    void
    operator()(const modificationExchange& rExchange) const
    {
      // Get the index, in the actual reactant species, of the mol in which
      // the modification exchange will happen.
      cpx::molSpec modMolSpecTr
	= rEnabling.translateMolSpec(rExchange.modMolSpec);

      // Determine what mol the modification will be swapped in.
      cpx::molSpec paradigmModMolSpec
	= rResultToResultParadigm.forward.molMap[modMolSpecTr];

      // Get the mol and make sure it's a mod-mol.
      cml::cptModMol* pModMol
	= dynamic_cast<cml::cptModMol*>(rResultParadigmMols[paradigmModMolSpec]);

      if(! pModMol)
	throw exchangedModXcpt(rResultParadigmMols[paradigmModMolSpec]);

      // Get the current modification state.
      cpx::modMolState exchangedMolState
	= pModMol->externState(rProductMolParams[paradigmModMolSpec]);

      // Substitute in the modification.
      exchangedMolState[rExchange.modSiteNdx] = rExchange.pReplacementMod;

      // Intern the new modMol state.
      const cpx::molParam nuMolParam
	= pModMol->internState(exchangedMolState);

      // Install the molParam into the vector of molParams for the primary
      // product species.
      rProductMolParams[paradigmModMolSpec] = nuMolParam;
    }
  };
  
  void
  omniRxnGen::
  respond(const fnd::featureStimulus<cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex> >& rStimulus)
  {
    // See if reaction generation is enabled.
    if(rCptUnit.getGenerateOk())
      {
	const cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex>& rNewContext
	  = rStimulus.getContext();
	
	// Construct the reaction and intern it for memory
	// management.
	cpt::globalReaction* pReaction
	  = new cpt::globalReaction(rCptUnit.getCompartmentGraph());
	pFamily->addEntry(pReaction);

	// Add the triggering complex as a reactant of multiplicity 1.
	pReaction->addReactant(rNewContext.getSpecies(),
			       1);

	// If an auxiliary reactant was specified, add it with multiplicity 1.
	if(pAdditionalReactant)
	  {
	    pReaction->addReactant(pAdditionalReactant,
				   1);
	  }

	// If an auxiliary product was specified, add it with multiplicity 1.
	if(pAdditionalProduct)
	  {
	    pReaction->addProduct(pAdditionalProduct,
				  1);
	  }

	// Set the rate of the reaction.
	pReaction->setRate(pExtrapolator->getRate(rNewContext));

	// Now construct the primary product species, the modified complex.
	// This occurs in two steps, making the structural complex and making
	// the parameters for the mols in the complex.

	// First, we have to construct its vector of mols, which is obtained
	// from the primary reactant's vector of mols by performing the
	// small-mol exchanges, if any.
	clx::cptPlex resultPlex(rNewContext.getPlexFamily().getParadigm());

	std::for_each(smallMolExchanges.begin(),
		      smallMolExchanges.end(),
		      exchangeSmallMols(rNewContext,
					resultPlex.mols));

	// Recognize the product complex, retaining the recognition map so
	// that we can reorder the states of the mols appropriately.
	cpx::plexIso resultToResultParadigm;
	clx::cptPlexFamily* pProductFamily
	  = rPlexUnit.recognize(resultPlex,
				resultToResultParadigm);

	// Now make the parameters for the mols in the product complex.

	// We start with a copy of the triggering complex's mol parameters.
	const std::vector<cpx::molParam>& rEnablingMolParams
	  = rNewContext.getMolParams();
	int paradigmMolNdx = rEnablingMolParams.size();
	std::vector<cpx::molParam> resultMolParams(paradigmMolNdx);
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
		      exchangeSmallMolParams(rNewContext,
					     resultToResultParadigm,
					     resultMolParams));

	// Perform the modification exchanges.
	std::for_each(modificationExchanges.begin(),
		      modificationExchanges.end(),
		      exchangeModifications(rNewContext,
					    resultToResultParadigm,
					    pProductFamily->getParadigm().mols,
					    resultMolParams));

	// Construct the primary product species.
	clx::cptPlexSpecies* pProductSpecies
	  = pProductFamily->getMember(resultMolParams);

	// Add the primary product species to the reaction, with multiplicity
	// 1.
	pReaction->addProduct(pProductSpecies,
			      1);

	// Continue reaction network generation at one lower depth.
	int notificationDepth = rStimulus.getNotificationDepth() - 1;
	if(0 <= notificationDepth)
	  pProductSpecies->ensureNotified(notificationDepth);

	// Finalize the reaction, creating and "hooking up" the
	// compartment reactions withe the compartment species.
	pReaction->finalizeCompartments(rCptApp.getPropensities());
      }
  }
}
