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

#ifndef DECOMPRXNGEN_H
#define DECOMPRXNGEN_H

/*! \file decompRxnGen.hh
  \ingroup decompGroup
  \brief Defines decomposition reaction generator. */

#include "mzr/reactionFamily.hh"
#include "plex/plexFeature.hh"
#include "plex/plexUnit.hh"
#include "dimer/decomposeExtrap.hh"

namespace dimer
{
  /*! \ingroup decompGroup
    \brief Reaction generator decompositions. */
  class decompRxnGen :
    public plx::bindingInContext::rxnGen
  {
    mzr::reactionFamily* pFamily;
    mzr::mzrUnit& rMzrUnit;
    plx::plexUnit& rPlexUnit;
    decomposeExtrapolator* pExtrap;
  
  public:
    // Note that this reaction generator memory manages the rate extrapolator.
    decompRxnGen(mzr::reactionFamily* pDecompFamily,
		 mzr::mzrUnit& refMzrUnit,
		 plx::plexUnit& refPlexUnit,
		 decomposeExtrapolator* pExtrapolator) :
      pFamily(pDecompFamily),
      rMzrUnit(refMzrUnit),
      rPlexUnit(refPlexUnit),
      pExtrap(pExtrapolator)
    {}
    
    ~decompRxnGen(void)
    {
      delete pExtrap;
    }

    void
    makeReactions(const plx::bindingInContext& rNewContext) const;
  };
}

#endif
