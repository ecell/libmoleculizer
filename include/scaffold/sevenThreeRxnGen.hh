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

#ifndef SEVENTHREERXNGEN_H
#define SEVENTHREERXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "plex/plexSpec.hh"
#include "plex/cxOmniParam.hh"
#include "scaffold/sevenThreeExtrap.hh"

namespace scaf
{
  class sevenThreeRxnGen :
    public plx::omniInContext::rxnGen
  {
    mzr::mzrUnit& rMzrUnit;
    
    // The kinase mol.
    bnd::modMol* pSte7;
    // Index of modification site at which AXP attaches to kinase.
    int atpModNdx;
    // Position of kinase in enabling omniplex.
    plx::plexMolSpec ste7Spec;
    // Index of modification site on Ste7 whose phosphorylation by Ste11
    // is a precondition for reaction.
    int activeModNdx;
    // Index of modification site on Ste7 whose phosphorylation by Fus3
    // inhibits reaction.
    int inhibModNdx;

    // The "substrate" mol.
    bnd::modMol* pFus3;
    // Position of substrate mol in enabling omniplex.
    plx::plexMolSpec fus3Spec;
    // Mask giving modification sites on Fus3 that can be phosphorylated
    // by Ste7.
    std::vector<bool> activityMask;

    const bnd::modification* pAtpBound;
    const bnd::modification* pAdpBound;
    const bnd::modification* pPhosphorylated;
    const bnd::modification* pNone;

    // This should be in the rxnGen base class template.
    mzr::reactionFamily* pFamily;

    // Reaction rate extrapolator.
    sevenThreeExtrapolator* pExtrap;

  public:

    sevenThreeRxnGen(mzr::mzrUnit& refMzrUnit,

		     bnd::modMol* pSte7ModMol,
		     int ste7AtpModSiteNdx,
		     plx::plexMolSpec ste7MolSpec,
		     int activePhosModSiteNdx,
		     int inhibPhosModSiteNdx,

		     bnd::modMol* pFus3ModMol,
		     plx::plexMolSpec fus3MolSpec,
		     const std::vector<bool>& rActivityMask,

		     const bnd::modification* pAtpBoundMod,
		     const bnd::modification* pAdpBoundMod,
		     const bnd::modification* pPhosphorylatedMod,
		     const bnd::modification* pNoneMod,

		     mzr::reactionFamily* pReactionFamily,
		     sevenThreeExtrapolator* pExtrapolator) :
      rMzrUnit(refMzrUnit),

      pSte7(pSte7ModMol),
      atpModNdx(ste7AtpModSiteNdx),
      ste7Spec(ste7MolSpec),
      activeModNdx(activePhosModSiteNdx),
      inhibModNdx(inhibPhosModSiteNdx),

      pFus3(pFus3ModMol),
      fus3Spec(fus3MolSpec),
      activityMask(rActivityMask),

      pAtpBound(pAtpBoundMod),
      pAdpBound(pAdpBoundMod),
      pPhosphorylated(pPhosphorylatedMod),
      pNone(pNoneMod),

      pFamily(pReactionFamily),
      pExtrap(pExtrapolator)
    {}

    ~sevenThreeRxnGen(void)
    {
      delete pExtrap;
    }

    void
    makeReactions(const plx::omniInContext& rContext) const;
  };
}

#endif // SEVENTHREERXNGEN_H
