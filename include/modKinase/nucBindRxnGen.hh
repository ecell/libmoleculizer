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

#ifndef NUCBINDRXNGEN_H
#define NUCBINDRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "stoch/stochSpecies.hh"
#include "modKinase/nucBindExtrap.hh"

namespace kinase
{
  // Basic nucleotide binding reaction.  This could have gone equally well
  // in the nucleotide exchange unit, but kinases and phosphatases are more
  // common, and using them would require basic nucleotide binding like this.
  //
  // This reaction generator makes both the binding and unbinding reactions.

  class nucBindRxnGen :
    public plx::molInContext::rxnGen
  {
    // The nucleotide-binding protein.
    bnd::modMol* pMol;

    // Index of the modification site representing nucleotide binding.
    int modSiteNdx;

    // The modifications denoting, respectively, nuclotide binding and nothing.
    const bnd::modification* pBound;
    const bnd::modification* pNone;

    // The nucleotide.
    stoch::stochSpecies* pNuc;

    // For the time being, this family will hold both the binding and unbinding
    // reactions.
    mzr::reactionFamily* pFamily;

    // When a reaction is sensitized to substrates, it is automaticaly
    // sensitized to the volume, which resides in mzr.
    mzr::mzrUnit& rMzrUnit;

    // Reaction rate extrapolators for the two kinds of reactions that
    // this reaction generator generates.
    nucBindExtrapolator* pNucBindExtrap;
    nucUnbindExtrapolator* pNucUnbindExtrap;
    
  public:
    nucBindRxnGen(bnd::modMol* pNucBindMol,
		  int nucModSiteNdx,
		  const bnd::modification* pNucBoundMod,
		  const bnd::modification* pNoneMod,
		  stoch::stochSpecies* pNucleotide,
		  mzr::reactionFamily* pReactionFamily,
		  mzr::mzrUnit& refMzrUnit,
		  nucBindExtrapolator* pNucBindExtrapolator,
		  nucUnbindExtrapolator* pNucUnbindExtrapolator) :
      pMol(pNucBindMol),
      modSiteNdx(nucModSiteNdx),
      pBound(pNucBoundMod),
      pNone(pNoneMod),
      pNuc(pNucleotide),
      pFamily(pReactionFamily),
      rMzrUnit(refMzrUnit),
      pNucBindExtrap(pNucBindExtrapolator),
      pNucUnbindExtrap(pNucUnbindExtrapolator)
    {}

    void
    makeReactions(const plx::molInContext& rContext) const;
  };
}

#endif // NUCBINDRXNGEN_H
