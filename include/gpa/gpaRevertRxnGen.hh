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

#ifndef GPAREVERTRXNGEN_H
#define GPAREVERTRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/modMol.hh"
#include "plex/plexFeature.hh"
#include "gpa/gpaRevertExtrap.hh"

namespace gpa
{
  class gpaRevertRxnGen :
    public plx::molInContext::rxnGen
  {
    mzr::reactionFamily* pFamily;

    // The gpa mol that will undergo reversion.
    bnd::modMol* pMol;

    // The index of the modification site where GDP/GTP bind.
    int gtpSiteNdx;

    // The enabling modification at the above site.
    const bnd::modification* pGtpBound;
    // The result modification at the above site.
    const bnd::modification* pGdpBound;

    // Product species.
    mzr::species* pPhosphate;

    // When a reaction is sensitized to its substrates, it has
    // to be sensitized to volume, which resides in mzrUnit.
    mzr::mzrUnit& rMzrUnit;

    // Rate extrapolator.
    gpaRevertExtrapolator* pExtrap;

  public:

    gpaRevertRxnGen(bnd::modMol* pGpaMol,
		    mzr::reactionFamily* pReactionFamily,
		    int gtpModSiteNdx,
		    const bnd::modification* pGdpBoundMod,
		    const bnd::modification* pGtpBoundMod,
		    mzr::species* pPhosphateSpecies,
		    mzr::mzrUnit& refMzrUnit,
		    gpaRevertExtrapolator* pExtrapolator) :
      pFamily(pReactionFamily),
      pMol(pGpaMol),
      gtpSiteNdx(gtpModSiteNdx),
      pGtpBound(pGtpBoundMod),
      pGdpBound(pGdpBoundMod),
      pPhosphate(pPhosphateSpecies),
      rMzrUnit(refMzrUnit),
      pExtrap(pExtrapolator)
    {}

    void
    makeReactions(const plx::molInContext& rContext) const;
  };
}

#endif // GPAREVERTRXNGEN_H
