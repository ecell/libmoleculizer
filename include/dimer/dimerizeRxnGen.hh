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

#ifndef DIMERIZERXNGEN_H
#define DIMERIZERXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mzr/binaryRxnGen.hh"
#include "plex/plexFeature.hh"
#include "plex/plexUnit.hh"
#include "dimer/dimerizeExtrap.hh"

namespace dimer
{
  class dimerUnit;
  
  /*! \file dimerizeRxnGen.hh
    \ingroup dimerizationGroup
    \brief Defines reaction generator for dimerizations. */

  /*! \ingroup dimerizationGroup
    \brief Reaction generators for dimerizations.

    There is typically one reaction generator for each substrate of the
    reactions generated. This class just simplifies the writing of both
    reaction generators for a binary reaction.  */
  class dimerizeRxnGenPair :
    public mzr::binaryRxnGenPair<plx::siteFeature, plx::siteFeature>
  {
    mzr::reactionFamily* pFamily;
    mzr::mzrUnit& rMzrUnit;
    plx::plexUnit& rPlexUnit;
    dimerizeExtrapolator* pExtrap;
    
  public:

    // Note that this reaction generator memory manages the rate extrapolator.
    dimerizeRxnGenPair(plx::siteFeature& rLeftSiteFeature,
		       plx::siteFeature& rRightSiteFeature,
		       mzr::reactionFamily* pDimerizeFamily,
		       mzr::mzrUnit& refMzrUnit,
		       plx::plexUnit& refPlexUnit,
		       dimerizeExtrapolator* pExtrapolator) :
      mzr::binaryRxnGenPair<plx::siteFeature, plx::siteFeature>
    (rLeftSiteFeature,
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
    makeBinaryReactions(const plx::siteInContext& rLeftContext,
			const plx::siteInContext& rRightContext) const;
  };
}

#endif
