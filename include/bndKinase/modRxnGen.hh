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

#ifndef MODRXNGEN_H
#define MODRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/smallMol.hh"
#include "plex/plexSpec.hh"
#include "plex/cxOmniParam.hh"
#include "bndKinase/modExtrap.hh"

namespace bndKinase
{
  class modRxnGen :
    public plx::omniInContext::rxnGen
  {
    // To intern reactions for memory management.
    mzr::mzrUnit& rMzrUnit;
    
    // To recognize the product complex.
    plx::plexUnit& rPlexUnit;

    // The substrate mod-mol.
    bnd::modMol* pSubstrate;

    // Index of the substrate mod-mol in the enabling omni.
    plx::plexMolSpec substrateSpec;

    // Index of the modification site on the substrate mod-mol.
    int substrateModSiteNdx;

    // Pointer to the modification to be installed.
    const bnd::modification* pInstalledMod;

    // Reaction family to receive the generated reactions.
    mzr::reactionFamily* pFamily;

    // Unary rate extrapolator.
    modExtrapolator* pExtrap;

  public:
    modRxnGen(mzr::mzrUnit& refMzrUnit,
	      plx::plexUnit& refPlexUnit,
	      bnd::modMol* pSubstrateModMol,
	      plx::plexMolSpec substrateMolSpec,
	      int modSiteNdx,
	      const bnd::modification* pInstalledMod,
	      mzr::reactionFamily* pReactionFamily,
	      modExtrapolator* pExtrapolator) :
      rMzrUnit(refMzrUnit),
      rPlexUnit(refPlexUnit),
      pSubstrate(pSubstrateModMol),
      substrateSpec(substrateMolSpec),
      substrateModSiteNdx(modSiteNdx),
      pInstalledMod(pInstalledMod),
      pFamily(pReactionFamily),
      pExtrap(pExtrapolator)
    {}

    ~modRxnGen(void)
    {
      delete pExtrap;
    }

    void
    makeReactions(const plx::omniInContext& rContext) const;
  };
}

#endif // MODRXNGEN_H
