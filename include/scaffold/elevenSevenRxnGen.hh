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

#ifndef ELEVENSEVENRXNGEN_H
#define ELEVENSEVENRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "plex/plexSpec.hh"
#include "plex/cxOmniParam.hh"
#include "scaffold/elevenSevenExtrap.hh"

namespace scaf
{
  class elevenSevenRxnGen :
    public plx::omniInContext::rxnGen
  {
    mzr::mzrUnit& rMzrUnit;
    
    // The kinase mol.
    bnd::modMol* pSte11;
    // Modification site where AXP attaches.
    int atpModSiteNdx;
    // Position in the enabling omniplex.
    plx::plexMolSpec ste11Spec;
    // The number of times Ste11 must be phosphorylated to be active as a
    // kinase of Ste7.  Right now, I'm counting all phosphorylation of Ste11,
    // not specified sites, though my notes indicate that N-terminal
    // phosphorylations are what count here.
    int elevenActiveCount;

    // The "substrate" mol.
    bnd::modMol* pSte7;
    // Position in the enabling omniplex.
    plx::plexMolSpec ste7Spec;
    // Site on Ste7 that can be phosphorylated by Ste11.  (This could easily
    // be generalized to a mask.)
    int targetModSiteNdx;

    const bnd::modification* pAtpBound;
    const bnd::modification* pAdpBound;
    const bnd::modification* pPhosphorylated;
    const bnd::modification* pNone;

    // This should be in the rxnGen base class template.
    mzr::reactionFamily* pFamily;

    // Unary rate extrapolator.
    elevenSevenExtrapolator* pExtrap;
    
  public:

    elevenSevenRxnGen(mzr::mzrUnit& refMzrUnit,

		      bnd::modMol* pSte11ModMol,
		      int ste11AtpModSiteNdx,
		      plx::plexMolSpec ste11MolSpec,
		      int elevenActivePhosCount,

		      bnd::modMol* pSte7ModMol,
		      plx::plexMolSpec ste7MolSpec,
		      int ste7TargetModSiteNdx,

		      const bnd::modification* pAtpBoundMod,
		      const bnd::modification* pAdpBoundMod,
		      const bnd::modification* pPhosphorylatedMod,
		      const bnd::modification* pNoneMod,

		      mzr::reactionFamily* pReactionFamily,
		      elevenSevenExtrapolator* pExtrapolator) :
      rMzrUnit(refMzrUnit),

      pSte11(pSte11ModMol),
      atpModSiteNdx(ste11AtpModSiteNdx),
      ste11Spec(ste11MolSpec),
      elevenActiveCount(elevenActivePhosCount),
      pSte7(pSte7ModMol),
      ste7Spec(ste7MolSpec),
      targetModSiteNdx(ste7TargetModSiteNdx),
  
      pAtpBound(pAtpBoundMod),
      pAdpBound(pAdpBoundMod),
      pPhosphorylated(pPhosphorylatedMod),
      pNone(pNoneMod),

      pFamily(pReactionFamily),
      pExtrap(pExtrapolator)
    {}

    ~elevenSevenRxnGen(void)
    {
      delete pExtrap;
    }

    void
    makeReactions(const plx::omniInContext& rContext) const;
  };
}

#endif
