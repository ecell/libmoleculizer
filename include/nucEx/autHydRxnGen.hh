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

#ifndef AUTHYDRXNGEN_H
#define AUTHYDRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "stoch/stochSpecies.hh"
#include "nucEx/autHydExtrap.hh"

namespace nucEx
{
  // Reaction generator for nucleotide auto-hydrolysis reaction, in which
  // a nucleotide-binding protein (Gpa1) hydrolyses bound, phosphorylated
  // nucleotide (GTP) releasing a phosphate.  It is much like a kinase
  // without a substrate.
  //
  // There is another reaction generator, hetHydRxnGen, for the
  // "hetero-hydrolysis" reaction, in which another enzyme (Sst2) assists in
  // performing the same reaction.  This reaction is also like some sort of
  // "broken GTP kinase," as though Sst2 were a substrate, but the phosphate
  // just doesn't get attached somehow.
  //
  // Alejandro says that there are no GTP-based kinases, though.

  class autHydRxnGen :
    public plx::molInContext::rxnGen
  {
    // The nucleotide-binding protein.
    bnd::modMol* pMol;

    // Index of the modification site representing nucleotide binding.
    int modSiteNdx;

    // The modifications denoting bound GTP/GDP. (Do ATP/ADP ever actually
    // appear in this context?
    const bnd::modification* pUnhydrolysedBound;
    const bnd::modification* pHydrolysedBound;

    // Phosphate is released as a side-product.
    stoch::stochSpecies* pPhosphate;

    mzr::reactionFamily* pFamily;

    // When a reaction is sensitized to substrates, it also must be
    // sensitized to volume, which resides in mzrUnit.
    mzr::mzrUnit& rMzrUnit;

    // Reaction rate extrapolator.
    autHydExtrapolator* pExtrap;
    
  public:
    autHydRxnGen(bnd::modMol* pAutHydMol,
		 int nucModSiteNdx,
		 const bnd::modification* pUnhydrolysedBoundMod,
		 const bnd::modification* pHydrolysedBoundMod,
		 stoch::stochSpecies* pPhosphateSpecies,
		 mzr::reactionFamily* pReactionFamily,
		 mzr::mzrUnit& refMzrUnit,
		 autHydExtrapolator* pExtrapolator) :
      pMol(pAutHydMol),
      modSiteNdx(nucModSiteNdx),
      pUnhydrolysedBound(pUnhydrolysedBoundMod),
      pHydrolysedBound(pHydrolysedBoundMod),
      pPhosphate(pPhosphateSpecies),
      pFamily(pReactionFamily),
      rMzrUnit(refMzrUnit),
      pExtrap(pExtrapolator)
    {}

    void
    makeReactions(const plx::molInContext& rContext) const;
  };
}

#endif // AUTHYDRXNGEN_H
