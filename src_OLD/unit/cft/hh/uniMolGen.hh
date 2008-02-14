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

#ifndef CFT_UNIMOLGEN_H
#define CFT_UNIMOLGEN_H

#include "cpx/cxMol.hh"
#include "cpt/cptUnit.hh"
#include "cml/cptModMol.hh"
#include "clx/clxUnit.hh"
#include "cft/uniMolExtrap.hh"

namespace cft
{
  class molModExchange
  {
  public:
    int modSiteNdx;
    const cpx::modification* pReplacementMod;

    molModExchange() :
      modSiteNdx(-1),
      pReplacementMod(0)
    {}

    molModExchange(int modSite,
		   const cpx::modification* pMod) :
      modSiteNdx(modSite),
      pReplacementMod(pMod)
    {}
  };

  class uniMolRxnGen :
    public fnd::rxnGen<cpx::cxMol<clx::cptPlexSpecies, clx::cptPlexFamily> >
  {
    // To add compartment reactions to the propensity distribution.
    cpt::cptApp& rCptApp;
    
    // To intern reactions for memory management.
    cpt::cptUnit& rCptUnit;

    // To recognize the new complex.
    clx::clxUnit& rClxUnit;

    // The mod-mol triggering reaction generation.
    cml::cptModMol* pEnablingMol;

    // Queries on the mol's state that must be passed to enable
    // reaction generation.
    cpx::andMolStateQueries* pMolQueries;

    // Exchanges of modifications.
    const std::vector<molModExchange> molModExchanges;

    // Additional massive reactant; null if there is no additional reactant.
    cpt::globalSpecies* pAdditionalReactant;

    // Additional product; null if there is no additional product.
    cpt::globalSpecies* pAdditionalProduct;

    // Reaction family to receive generated reactions.
    utl::autoVector<cpt::globalReaction>* pFamily;

    // Rate extrapolator, which sometimes is unary and sometimes is binary.
    // This omniRxnGen is responsible for memory management of its
    // extrapolator.
    const uniMolExtrapolator* pExtrapolator;

  public:
    uniMolRxnGen(cpt::cptApp& refCptApp,
		 cpt::cptUnit& refCptUnit,
		 clx::clxUnit& refClxUnit,
		 cml::cptModMol* pEnablingModMol,
		 cpx::andMolStateQueries* pAndMolStateQueries,
		 const std::vector<molModExchange>& rModExchanges,
		 cpt::globalSpecies* pAuxReactant,
		 cpt::globalSpecies* pAuxProduct,
		 utl::autoVector<cpt::globalReaction>* pReactionFamily,
		 const uniMolExtrapolator* pUniMolExtrapolator) :
      rCptApp(refCptApp),
      rCptUnit(refCptUnit),
      rClxUnit(refClxUnit),
      pEnablingMol(pEnablingModMol),
      pMolQueries(pAndMolStateQueries),
      molModExchanges(rModExchanges),
      pAdditionalReactant(pAuxReactant),
      pAdditionalProduct(pAuxProduct),
      pFamily(pReactionFamily),
      pExtrapolator(pUniMolExtrapolator)
    {}

    ~uniMolRxnGen(void)
    {
      delete pExtrapolator;
    }

    void
    respond(const fnd::featureStimulus<cpx::cxMol<clx::cptPlexSpecies, clx::cptPlexFamily> >& rStimulus);
  };
}

#endif // CFT_UNIMOLGEN_H
