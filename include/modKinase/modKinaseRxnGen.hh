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

#ifndef MODKINASERXNGEN_H
#define MODKINASERXNGEN_H

#include "mol/modMol.hh"
#include "mzr/reactionFamily.hh"
#include "mzr/binaryRxnGen.hh"
#include "plex/plexFeature.hh"
#include "modKinase/modKinaseExtrap.hh"

namespace kinase
{
  class modKinaseRxnGenPair :
    public mzr::binaryRxnGenPair<plx::molFeature, plx::molFeature>
  {
    bnd::modMol* pKinase;
    bnd::modMol* pSubstrate;
    mzr::reactionFamily* pFamily;

    // The index of the ATP binding site on the kinase.
    int atpSiteNdx;

    // This mask tells which sites on the substrate can be phosphorylated
    // by the kinase.
    std::vector<bool> activityMask;

    // As a temporary stopgap, I'll look these up by name in the
    // constructor.
    const bnd::modification* pPhosphorylated;
    const bnd::modification* pATPbound;
    const bnd::modification* pADPbound;
    const bnd::modification* pNone;

    // When a reaction is sensitized to its substrates, it is automatically
    // sensitized to volume, which resides in mzr.
    mzr::mzrUnit& rMzrUnit;

    // Rate extrapolator.
    modKinaseExtrapolator* pExtrap;

  public:

    modKinaseRxnGenPair(bnd::modMol* pKinaseMol,
			int kinaseAtpBindingSiteIndex,
			bnd::modMol* pSubstrateMol,
			const bnd::modification* pPhosphorylatedMod,
			const bnd::modification* pATPboundMod,
			const bnd::modification* pADPboundMod,
			const bnd::modification* pNoneMod,
			const std::vector<bool>& rActivityMask,
			mzr::reactionFamily* pReactionFamily,
			mzr::mzrUnit& refMzrUnit,
			modKinaseExtrapolator* pExtrapolator) :
      mzr::binaryRxnGenPair<plx::molFeature, plx::molFeature>(*pKinaseMol,
					       *pSubstrateMol),
      pKinase(pKinaseMol),
      pSubstrate(pSubstrateMol),
      pFamily(pReactionFamily),
      atpSiteNdx(kinaseAtpBindingSiteIndex),
      activityMask(rActivityMask),
      pPhosphorylated(pPhosphorylatedMod),
      pATPbound(pATPboundMod),
      pADPbound(pADPboundMod),
      pNone(pNoneMod),
      rMzrUnit(refMzrUnit),
      pExtrap(pExtrapolator)
    {}

    ~modKinaseRxnGenPair(void)
    {
      delete pExtrap;
    }

    void
    makeBinaryReactions(const plx::molInContext& rKinaseContext,
			const plx::molInContext& rSubstrateContext) const;
  };
}

#endif // MODKINASERXNGEN_H
