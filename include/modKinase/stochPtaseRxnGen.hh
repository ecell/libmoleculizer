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

#ifndef STOCHPTASERXNGEN_H
#define STOCHPTASERXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "stoch/stochSpecies.hh"
#include "modKinase/stochPtaseExtrap.hh"

namespace kinase
{
  class stochPtaseRxnGen :
    public plx::molInContext::rxnGen
  {
    // The substrate.
    //
    // This mol will be acted upon by the phosphatase regardless of its
    // molecular context.
    bnd::modMol* pMol;

    // The phosphatase.
    stoch::stochSpecies* pPtase;

    // Phosphate, one of the reaction products.
    stoch::stochSpecies* pPhosphate;

    // Modification representing the phosphorylated state.
    const bnd::modification* pPhosMod;
    // Modification representing no modification, which is left
    // after dephosphorylation.
    const bnd::modification* pNoneMod;

    // Mask identifying the modification sites that can be dephosphorylated by
    // the phosphatase.
    std::vector<bool> phosMask;

    mzr::reactionFamily* pFamily;

    // When a reaction is sensitized to it substrates, it is automatically
    // sensitized to volume, which resides in mzr.
    mzr::mzrUnit& rMzrUnit;

    stochPtaseExtrapolator* pExtrap;
    
  public:

    stochPtaseRxnGen(bnd::modMol* pSubstrateModMol,
		     stoch::stochSpecies* pStochPtase,
		     stoch::stochSpecies* pStochPhosphate,
		     const bnd::modification* pPhosModification,
		     const bnd::modification* pNoneModification,
		     const std::vector<bool>& rModSiteMask,
		     mzr::reactionFamily* pReactionFamily,
		     mzr::mzrUnit& refMzrUnit,
		     stochPtaseExtrapolator* pExtrapolator) :
      pMol(pSubstrateModMol),
      pPtase(pStochPtase),
      pPhosphate(pStochPhosphate),
      pPhosMod(pPhosModification),
      pNoneMod(pNoneModification),
      phosMask(rModSiteMask),
      pFamily(pReactionFamily),
      rMzrUnit(refMzrUnit),
      pExtrap(pExtrapolator)
    {}

    void
    makeReactions(const plx::molInContext& rContext) const;
  };
}

#endif // STOCHPTASERXNGEN_H
