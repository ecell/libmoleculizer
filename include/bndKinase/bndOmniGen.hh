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

#ifndef BNDOMNIGEN_H
#define BNDOMNIGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/smallMol.hh"
#include "plex/plexSpec.hh"
#include "plex/plexFeature.hh"
#include "plex/cxOmniParam.hh"
#include "bndKinase/bndOmniExtrap.hh"

namespace bndKinase
{
  class smallMolExchange
  {
  public:
    plx::plexMolSpec exchangedMolSpec;
    bnd::smallMol* pReplacementMol;

    smallMolExchange(void) :
      exchangedMolSpec(-1),
      pReplacementMol(0)
    {}
    
    smallMolExchange(plx::plexMolSpec xchngMolSpec,
		     bnd::smallMol* pRplcmntMol) :
      exchangedMolSpec(xchngMolSpec),
      pReplacementMol(pRplcmntMol)
    {}
  };

  class modificationExchange
  {
  public:
    plx::plexMolSpec modMolSpec;
    int modSiteNdx;
    const bnd::modification* pReplacementMod;

    modificationExchange() :
      modMolSpec(-1),
      modSiteNdx(-1),
      pReplacementMod(0)
    {}

    modificationExchange(plx::plexMolSpec molSpec,
			 int modSite,
			 const bnd::modification* pMod) :
      modMolSpec(molSpec),
      modSiteNdx(modSite),
      pReplacementMod(pMod)
    {}
  };

  class bndOmniRxnGen :
    public plx::omniInContext::rxnGen
  {
    // To intern reactions for memory management.
    mzr::mzrUnit& rMzrUnit;

    // To recognize the new complex.
    plx::plexUnit& rPlexUnit;

    // Exchanges of small-mol components.
    const std::vector<smallMolExchange> smallMolExchanges;

    // Exchanges of modifications.
    const std::vector<modificationExchange> modificationExchanges;

    // Additional reactant; null if there is no additional reactant.
    mzr::species* pAdditionalReactant;

    // Additional product; null if there is no additional product.
    mzr::species* pAdditionalProduct;

    // Reaction family to receive generated reactions.
    mzr::reactionFamily* pFamily;

    // Rate extrapolator, which sometimes is unary and sometimes is binary.
    // This bndOmniRxnGen is responsible for memory management of its
    // extrapolator.
    const bndOmniExtrapolator* pExtrapolator;

  public:
    bndOmniRxnGen(mzr::mzrUnit& refMzrUnit,
		  plx::plexUnit& refPlexUnit,
		  const std::vector<smallMolExchange>& rSMExchanges,
		  const std::vector<modificationExchange>& rModExchanges,
		  mzr::species* pAuxReactant,
		  mzr::species* pAuxProduct,
		  mzr::reactionFamily* pReactionFamily,
		  const bndOmniExtrapolator* pBndOmniExtrapolator) :
      rMzrUnit(refMzrUnit),
      rPlexUnit(refPlexUnit),
      smallMolExchanges(rSMExchanges),
      modificationExchanges(rModExchanges),
      pAdditionalReactant(pAuxReactant),
      pAdditionalProduct(pAuxProduct),
      pFamily(pReactionFamily),
      pExtrapolator(pBndOmniExtrapolator)
    {}

    ~bndOmniRxnGen(void)
    {
      delete pExtrapolator;
    }

    void
    makeReactions(const plx::omniInContext& rContext) const;
  };
}

#endif //  BNDOMNIGEN_H
