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

#ifndef GPAEXCHANGERXNGEN_H
#define GPAEXCHANGERXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "plex/plexSpec.hh"
#include "stoch/stochSpecies.hh"
#include "gpa/gpaExchangeExtrap.hh"

namespace gpa
{
  class gpaExchangeRxnGen :
    public plx::omniInContext::rxnGen
  {
    // This should probably be added to the various rxnGen base classes.
    mzr::reactionFamily* pFamily;

    // The index of the gpa mol in the enabling receptor complex.
    plx::plexMolSpec gpaSpec;

    // The gpa mol itself.
    bnd::modMol* pMol;

    // Index of the modification site where GDP/GTP bind.
    int gtpSiteNdx;

    // The enabling modification at the above site.
    const bnd::modification* pGdpBound;
    // The result modificaiton at the above site.
    const bnd::modification* pGtpBound;

    // A product species.  Actually, we only need this as an abstract
    // species.
    stoch::stochSpecies* pGDP;
    // A substrate species.  We need this one to be massive,
    // in order to convert the weight-independent reaction constant
    // to the actual reaction rates for generated reactions.
    stoch::stochSpecies* pGTP;

    // This is for memory management of reactions; used in
    // reactionFamily::addReaction.
    mzr::mzrUnit& rMzrUnit;

    // Rate extrapolator.
    gpaExchangeExtrapolator* pExtrap;

  public:

    gpaExchangeRxnGen(bnd::modMol* pGpaMol,
		      plx::plexMolSpec gpaMolSpec,
		      int gtpModSiteNdx,
		      const bnd::modification* pGdpBoundMod,
		      const bnd::modification* pGtpBoundMod,
		      stoch::stochSpecies* pGdpSpecies,
		      stoch::stochSpecies* pGtpSpecies,
		      mzr::reactionFamily* pReactionFamily,
		      mzr::mzrUnit& refMzrUnit,
		      gpaExchangeExtrapolator* pExtrapolator) :
      pFamily(pReactionFamily),
      gpaSpec(gpaMolSpec),
      pMol(pGpaMol),
      gtpSiteNdx(gtpModSiteNdx),
      pGdpBound(pGdpBoundMod),
      pGtpBound(pGtpBoundMod),
      pGDP(pGdpSpecies),
      pGTP(pGtpSpecies),
      rMzrUnit(refMzrUnit),
      pExtrap(pExtrapolator)
    {}

    ~gpaExchangeRxnGen(void)
    {
      delete pExtrap;
    }

    void
    makeReactions(const plx::omniInContext& rContext) const;
  };
}

#endif
