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

#ifndef BNDKINASERXNGEN_H
#define BNDKINASERXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/smallMol.hh"
#include "plex/plexSpec.hh"
#include "plex/cxOmniParam.hh"
#include "bndKinase/bndKinaseExtrap.hh"

namespace bndKinase
{
  class bndKinaseRxnGen :
    public plx::omniInContext::rxnGen
  {
    // To intern reactions for memory management.
    mzr::mzrUnit& rMzrUnit;

    // To recognize the new complex.
    plx::plexUnit& rPlexUnit;
    
    // The substrate mod-mol.
    bnd::modMol* pSubstrate;

    // Index of the substrate mol.
    plx::plexMolSpec substrateSpec;

    // Index of the ATP mol in the omni.
    plx::plexMolSpec atpSpec;

    // Pointer to the ADP mol, with which to replace the ATP mol.
    bnd::smallMol* pADPmol;

    // Index of the phosphorylation modification site on the substrate.
    int phosSiteNdx;

    // Pointer to the phosphorylated modification.
    const bnd::modification* pPhosphorylated;

    // Reaction family to receive generated reactions.
    mzr::reactionFamily* pFamily;

    // Unary rate extrapolator.
    bndKinaseExtrapolator* pExtrap;

  public:
    bndKinaseRxnGen(mzr::mzrUnit& refMzrUnit,
		    plx::plexUnit& refPlexUnit,
		    bnd::modMol* pSubstrateModMol,
		    plx::plexMolSpec substrateMolSpec,
		    plx::plexMolSpec atpMolSpec,
		    bnd::smallMol* pADPsmallMol,
		    int phosModSiteNdx,
		    const bnd::modification* pPhosphorylatedMod,
		    mzr::reactionFamily* pReactionFamily,
		    bndKinaseExtrapolator* pExtrapolator) :
      rMzrUnit(refMzrUnit),
      rPlexUnit(refPlexUnit),
      pSubstrate(pSubstrateModMol),
      substrateSpec(substrateMolSpec),
      atpSpec(atpMolSpec),
      pADPmol(pADPsmallMol),
      phosSiteNdx(phosModSiteNdx),
      pPhosphorylated(pPhosphorylatedMod),
      pFamily(pReactionFamily),
      pExtrap(pExtrapolator)
    {}

    ~bndKinaseRxnGen(void)
    {
      delete pExtrap;
    }

    void
    makeReactions(const plx::omniInContext& rContext) const;
  };
}

#endif // BNDKINASERXNGEN_H
