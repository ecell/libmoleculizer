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

#ifndef CFT_OMNIGEN_H
#define CFT_OMNIGEN_H

#include "cpx/cxOmni.hh"
#include "cml/cptSmallMol.hh"
#include "cft/omniExtrap.hh"

namespace cft
{
  class smallMolExchange
  {
  public:
    cpx::molSpec exchangedMolSpec;
    cml::cptSmallMol* pReplacementMol;

    smallMolExchange(void) :
      exchangedMolSpec(-1),
      pReplacementMol(0)
    {}
    
    smallMolExchange(cpx::molSpec xchngMolSpec,
		     cml::cptSmallMol* pRplcmntMol) :
      exchangedMolSpec(xchngMolSpec),
      pReplacementMol(pRplcmntMol)
    {}
  };

  class modificationExchange
  {
  public:
    cpx::molSpec modMolSpec;
    int modSiteNdx;
    const cpx::modification* pReplacementMod;

    modificationExchange() :
      modMolSpec(-1),
      modSiteNdx(-1),
      pReplacementMod(0)
    {}

    modificationExchange(cpx::molSpec molSpec,
			 int modSite,
			 const cpx::modification* pMod) :
      modMolSpec(molSpec),
      modSiteNdx(modSite),
      pReplacementMod(pMod)
    {}
  };

  class omniRxnGen :
    public fnd::rxnGen<cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex> >
  {
    // To add compartment reactions to the propensity distribution.
    cpt::cptApp& rCptApp;
    
    // To intern reactions for memory management.
    cpt::cptUnit& rCptUnit;

    // To recognize the new complex.
    clx::clxUnit& rPlexUnit;

    // Exchanges of small-mol components.
    const std::vector<smallMolExchange> smallMolExchanges;

    // Exchanges of modifications.
    const std::vector<modificationExchange> modificationExchanges;

    // Additional reactant; null if there is no additional reactant.
    cpt::globalSpecies* pAdditionalReactant;

    // Additional product; null if there is no additional product.
    cpt::globalSpecies* pAdditionalProduct;

    // Reaction family to receive generated reactions.
    utl::autoVector<cpt::globalReaction>* pFamily;

    // Rate extrapolator, which sometimes is unary and sometimes is binary.
    // This omniRxnGen is responsible for memory management of its
    // extrapolator.
    const omniExtrapolator* pExtrapolator;

  public:
    omniRxnGen(cpt::cptApp& refCptApp,
	       cpt::cptUnit& refCptUnit,
	       clx::clxUnit& refPlexUnit,
	       const std::vector<smallMolExchange>& rSMExchanges,
	       const std::vector<modificationExchange>& rModExchanges,
	       cpt::globalSpecies* pAuxReactant,
	       cpt::globalSpecies* pAuxProduct,
	       utl::autoVector<cpt::globalReaction>* pReactionFamily,
	       const omniExtrapolator* pOmniExtrapolator) :
      rCptApp(refCptApp),
      rCptUnit(refCptUnit),
      rPlexUnit(refPlexUnit),
      smallMolExchanges(rSMExchanges),
      modificationExchanges(rModExchanges),
      pAdditionalReactant(pAuxReactant),
      pAdditionalProduct(pAuxProduct),
      pFamily(pReactionFamily),
      pExtrapolator(pOmniExtrapolator)
    {}

    ~omniRxnGen(void)
    {
      delete pExtrapolator;
    }

    void
    respond(const fnd::featureStimulus<cpx::cxOmni<cml::cptMol, clx::cptPlexSpecies, clx::cptPlexFamily, clx::cptOmniPlex> >& rStimulus);
  };
}

#endif //  CFT_OMNIGEN_H
