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

#ifndef HETHYDRXNGEN_H
#define HETHYDRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "stoch/stochSpecies.hh"
#include "nucEx/hetHydExtrap.hh"

namespace nucEx
{
  // Reaction generator for a nuclotide hydrolysis reaction that is "assisted"
  // by an enzyme (Sst2) other than the one undergoing the hydrolysis reaction
  // (Gpa1).
  //
  // (The reaction generator autHydRxnGen was designed to model the
  // "auto-reversion" reaction of Gpa1.)

  class hetHydRxnGen :
    public plx::omniInContext::rxnGen
  {
    // The nucleotide-binding protein.
    bnd::modMol* pMol;

    // Index of the nucleotide-binding protein in the "enabling" subcomplex,
    // consisting of the n-b protein and the "assistant" protien.
    plx::plexMolSpec nbMolSpec;

    // Index of the site representing nucleotide binding.
    int modSiteNdx;

    // The modifications denoting bound GTP/GDP. (Do ATP/ADP ever actually
    // appear in this context?
    const bnd::modification* pUnhydrolysedBound;
    const bnd::modification* pHydrolysedBound;

    // Phosphate is released as a side-product.
    stoch::stochSpecies* pPhosphate;

    mzr::reactionFamily* pFamily;

    // Reations must be sensitized to volume, which resides in mzrUnit.
    mzr::mzrUnit& rMzrUnit;
    
    // Unary reaction rate extrapolator. (The complex containing the
    // "assistant" mol and the nucleotide-binding mol is the only substrate.)
    hetHydExtrapolator* pExtrap;

  public:
    hetHydRxnGen(bnd::modMol* pTargetMol,
		 int nucBindingMolSpec,
		 int nucModSiteNdx,
		 const bnd::modification* pUnhydrolysedBoundMod,
		 const bnd::modification* pHydrolysedBoundMod,
		 stoch::stochSpecies* pPhosphateSpecies,
		 mzr::reactionFamily* pReactionFamily,
		 mzr::mzrUnit& refMzrUnit,
		 hetHydExtrapolator* pExtrapolator) :
      pMol(pTargetMol),
      nbMolSpec(nucBindingMolSpec),
      modSiteNdx(nucModSiteNdx),
      pUnhydrolysedBound(pUnhydrolysedBoundMod),
      pHydrolysedBound(pHydrolysedBoundMod),
      pPhosphate(pPhosphateSpecies),
      pFamily(pReactionFamily),
      rMzrUnit(refMzrUnit),
      pExtrap(pExtrapolator)
    {}

    ~hetHydRxnGen(void)
    {
      delete pExtrap;
    }

    void
    makeReactions(const plx::omniInContext& rContext) const;
  };
}

#endif // HETHYDRXNGEN_H
