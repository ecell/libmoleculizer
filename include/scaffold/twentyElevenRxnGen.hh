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

#ifndef TWENTYELEVENRXNGEN_H
#define TWENTYELEVENRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "mol/molUnit.hh"
#include "plex/plexFeature.hh"
#include "plex/plexSpec.hh"
#include "plex/cxOmniParam.hh"
#include "scaffold/twentyElevenExtrap.hh"

namespace scaf
{
  // This is almost like an ordinary kinase generator, except that reactions are
  // only generated when both the substrate mol and the kinase mol are in a
  // particular enabling sub-complex.  Thus, this reaction generator listens for
  // appearances of complexes containing the enabling sub-complex.
  class twentyElevenRxnGen :
    public plx::omniInContext::rxnGen
  {
    mzr::mzrUnit& rMzrUnit;
    
    // The kinase and its info.
    bnd::modMol* pSte20;
    int atpModSiteNdx;
    // The index in the enabling omniplex.
    plx::plexMolSpec ste20Spec;

    // The substrate and its info.
    bnd::modMol* pSte11;
    // This will always be necessary just to distinguish
    // phosphorylation sites from other modification sites.
    std::vector<bool> activityMask;
    // The index in the enabling omniplex.
    plx::plexMolSpec ste11Spec;

    const bnd::modification* pAtpBound;
    const bnd::modification* pAdpBound;
    const bnd::modification* pPhosphorylated;
    const bnd::modification* pNone;

    // This should probably be in rxnGen base class template.
    mzr::reactionFamily* pFamily;

    // Reaction rate extrapolator.
    twentyElevenExtrapolator* pExtrap;
    
  public:

    twentyElevenRxnGen(mzr::mzrUnit& refMzrUnit,

		       bnd::modMol* pSte20ModMol,
		       int ste20AtpBindingSiteNdx,
		       plx::plexMolSpec ste20MolSpec,

		       bnd::modMol* pSte11ModMol,
		       const std::vector<bool>& rActivityMask,
		       plx::plexMolSpec ste11MolSpec,

		       const bnd::modification* pAtpBoundMod,
		       const bnd::modification* pAdpBoundMod,
		       const bnd::modification* pPhosphorylatedMod,
		       const bnd::modification* pNoneMod,

		       mzr::reactionFamily* pReactionFamily,
		       twentyElevenExtrapolator* pExtrapolator) :
      rMzrUnit(refMzrUnit),
      
      pSte20(pSte20ModMol),
      atpModSiteNdx(ste20AtpBindingSiteNdx),
      ste20Spec(ste20MolSpec),

      pSte11(pSte11ModMol),
      activityMask(rActivityMask),
      ste11Spec(ste11MolSpec),

      pAtpBound(pAtpBoundMod),
      pAdpBound(pAdpBoundMod),
      pPhosphorylated(pPhosphorylatedMod),
      pNone(pNoneMod),

      pFamily(pReactionFamily),
      pExtrap(pExtrapolator)
    {}

    ~twentyElevenRxnGen(void)
    {
      delete pExtrap;
    }

    void
    makeReactions(const plx::omniInContext& rContext) const;
  };
}

#endif
