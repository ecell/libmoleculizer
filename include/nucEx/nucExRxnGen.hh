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

#ifndef NUCEXRXNGEN_H
#define NUCEXRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "stoch/stochSpecies.hh"
#include "nucEx/nucExExtrap.hh"

namespace nucEx
{
  // This reaction generator produces both nucleotide-binding and unbinding
  // reactions, with two sets of reaction rates, depending on whether the
  // binding mol is appropriately situated in the "enabling complex."
  // (Otherwise, this should be very similar to the basic nucleotide-binding/
  // unbinding reaction generator.)
  //
  // For g-protein coupled receptor, there will be two of these reaction
  // generators, one for GTP and one for GDP.
  class nucExRxnGen :
    public plx::molInContext::rxnGen
  {
    // The nucleotide-binding protein.
    bnd::modMol* pMol;

    // Index of the modification site being used to represent nucleotide
    // binding.
    int modSiteNdx;

    // This modification enables the binding reaction and is the result
    // of the unbinding reaction.
    const bnd::modification* pNoneMod;

    // This modification enables the unbinding reaction and is the result
    // of the binding reaction.
    const bnd::modification* pBoundMod;

    // The plexFamily of the enabling sub-complex, and the index of the
    // target mol in that plexFamily's paradigm.
    plx::plexFamily* pEnablingFamily;
    plx::plexMolSpec enablingSpec;

    // The nucleotide.
    stoch::stochSpecies* pNucSpecies;

    // The weight-independent version of the on-rates for the enabled
    // and non-enabled situations.
//     const double enabledOnConst;
//     const double plainOnConst;

    // The off-rates for the enabled and non-enabled situations.
//     const double enabledOffRate;
//     const double plainOffRate;

    // For now, put all the reactions into one reaction family.
    mzr::reactionFamily* pFamily;

    // When a reaction is sensitized to its substrates, it must
    // be sensitized as well to volume, which resides in mzrUnit.
    mzr::mzrUnit& rMzrUnit;

    // Reaction rate extrapolators for the four different kinds
    // of reactions produced by this reaction generator.
    nucExBindExtrapolator* pEnabledBindExtrap;
    nucExUnbindExtrapolator* pEnabledUnbindExtrap;
    
    nucExBindExtrapolator* pPlainBindExtrap;
    nucExUnbindExtrapolator* pPlainUnbindExtrap;
    
  public:
    nucExRxnGen(bnd::modMol* pTargetMol,
		int nucModSiteNdx,
		const bnd::modification* pNoneModification,
		const bnd::modification* pBoundModification,
		plx::plexFamily* pEnablingSubcomplex,
		plx::plexMolSpec enablingMolSpec,
		stoch::stochSpecies* pNucleotideSpecies,
		mzr::reactionFamily* pReactionFamily,
		mzr::mzrUnit& refMzrUnit,
		nucExBindExtrapolator* pEnabledBindExtrapolator,
		nucExUnbindExtrapolator* pEnabledUnbindExtrapolator,
		nucExBindExtrapolator* pPlainBindExtrapolator,
		nucExUnbindExtrapolator* pPlainUnbindExtrapolator) :
      pMol(pTargetMol),
      modSiteNdx(nucModSiteNdx),
      pNoneMod(pNoneModification),
      pBoundMod(pBoundModification),
      pEnablingFamily(pEnablingSubcomplex),
      enablingSpec(enablingMolSpec),
      pNucSpecies(pNucleotideSpecies),
      pFamily(pReactionFamily),
      rMzrUnit(refMzrUnit),
      pEnabledBindExtrap(pEnabledBindExtrapolator),
      pEnabledUnbindExtrap(pEnabledUnbindExtrapolator),
      pPlainBindExtrap(pPlainBindExtrapolator),
      pPlainUnbindExtrap(pPlainUnbindExtrapolator)
    {}

    ~nucExRxnGen(void)
    {
      delete pEnabledBindExtrap;
      delete pEnabledUnbindExtrap;
      delete pPlainBindExtrap;
      delete pPlainUnbindExtrap;
    }

    void
    makeReactions(const plx::molInContext& rContext) const;
  };
}

#endif // NUCEXRXNGEN_H
