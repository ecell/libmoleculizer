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

#ifndef DIMER_DIMERIZERXNGEN_H
#define DIMER_DIMERIZERXNGEN_H

#include "fnd/binaryRxnGen.hh"
#include "mol/siteFeature.hh"
#include "plex/plexUnit.hh"
#include "dimer/dimerizeExtrap.hh"

namespace dimer
{
  class dimerUnit;
  
  class dimerizeRxnGenPair :
    public fnd::binaryRxnGenPair<bnd::siteFeature, bnd::siteFeature>
  {
    utl::autoVector<mzr::mzrReaction>* pFamily;
    mzr::mzrUnit& rMzrUnit;
    plx::plexUnit& rPlexUnit;
    dimerizeExtrapolator* pExtrap;
    
  public:

    // Note that this reaction generator memory manages the rate extrapolator.
    dimerizeRxnGenPair(bnd::siteFeature& rLeftSiteFeature,
		       bnd::siteFeature& rRightSiteFeature,
		       utl::autoVector<mzr::mzrReaction>* pDimerizeFamily,
		       mzr::mzrUnit& refMzrUnit,
		       plx::plexUnit& refPlexUnit,
		       dimerizeExtrapolator* pExtrapolator) :
      fnd::binaryRxnGenPair<bnd::siteFeature, bnd::siteFeature>(rLeftSiteFeature,
								rRightSiteFeature),
      pFamily(pDimerizeFamily),
      rMzrUnit(refMzrUnit),
      rPlexUnit(refPlexUnit),
      pExtrap(pExtrapolator)
    {}

    ~dimerizeRxnGenPair(void)
    {
      delete pExtrap;
    }

    void
    makeBinaryReactions
    (const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rLeftContext,
     const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rRightContext,
     int generateDepth) const;
  };
}

#endif // DIMER_DIMERIZERXNGEN_H
