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

#ifndef THREESEVENRXNGEN_H
#define THREESEVENRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "plex/plexSpec.hh"
#include "plex/cxOmniParam.hh"
#include "scaffold/threeSevenExtrap.hh"

namespace scaf
{
  class threeSevenRxnGen :
    public plx::omniInContext::rxnGen
  {
    mzr::mzrUnit& rMzrUnit;
    
    // The kinase mol (Fus3).
    bnd::modMol* pFus3;
    // Modification site at which AXP attaches.
    int atpModNdx;
    // Position of kinase in enabling omniplex.
    plx::plexMolSpec fus3Spec;
    // Minimum number of times the kinase (Fus3) must be phosphorylated to be
    // active.
    int minFus3PhosCount;

    // The "substrate" mol (Ste7).
    bnd::modMol* pSte7;
    // Position in the enabling omniplex.
    plx::plexMolSpec ste7Spec;
    // Mask indicating which modification sites on Ste7 can be phosphorylated
    // by Fus3.
    std::vector<bool> activityMask;

    const bnd::modification* pAtpBound;
    const bnd::modification* pAdpBound;
    const bnd::modification* pPhosphorylated;
    const bnd::modification* pNone;

    // This should be in the rxnGen base class template?
    mzr::reactionFamily* pFamily;

    // The (unary) reaction rate.
    double rate;

    // Reaction rate extrapolator.
    threeSevenExtrapolator* pExtrap;

  public:

    threeSevenRxnGen(mzr::mzrUnit& refMzrUnit,

		     bnd::modMol* pFus3ModMol,
		     int fus3AtpModSiteNdx,
		     plx::plexMolSpec fus3MolSpec,
		     int minFus3ActivePhosCount,

		     bnd::modMol* pSte7ModMol,
		     plx::plexMolSpec ste7MolSpec,
		     std::vector<bool>& rActivityMask,

		     const bnd::modification* pAtpBoundMod,
		     const bnd::modification* pAdpBoundMod,
		     const bnd::modification* pPhosphorylatedMod,
		     const bnd::modification* pNoneMod,

		     mzr::reactionFamily* pReactionFamily,
		     threeSevenExtrapolator* pExtrapolator) :
      rMzrUnit(refMzrUnit),

      pFus3(pFus3ModMol),
      atpModNdx(fus3AtpModSiteNdx),
      fus3Spec(fus3MolSpec),
      minFus3PhosCount(minFus3ActivePhosCount),

      pSte7(pSte7ModMol),
      ste7Spec(ste7MolSpec),
      activityMask(rActivityMask),

      pAtpBound(pAtpBoundMod),
      pAdpBound(pAdpBoundMod),
      pPhosphorylated(pPhosphorylatedMod),
      pNone(pNoneMod),

      pFamily(pReactionFamily),
      pExtrap(pExtrapolator)
    {}

    ~threeSevenRxnGen(void)
    {
      delete pExtrap;
    }

    void
    makeReactions(const plx::omniInContext& rContext) const;
  };
}

#endif // THREESEVENRXNGEN_H
