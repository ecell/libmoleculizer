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

#ifndef TWENTYELEVENFAM_H
#define TWENTYELEVENFAM_H

#include "mzr/reactionFamily.hh"
#include "scaffold/twentyElevenRxnGen.hh"

namespace scaf
{
  class twentyElevenFam : public mzr::reactionFamily
  {
    twentyElevenRxnGen generator;
  
  public:

    twentyElevenFam(mzr::mzrUnit& rMzrUnit,

		    bnd::modMol* pSte20ModMol,
		    int ste20AtpBindingSiteNdx,
		    plx::plexMolSpec ste20MolSpec,

		    bnd::modMol* pSte11ModMol,
		    const std::vector<bool>& rActivityMask,
		    plx::plexMolSpec ste11MolSpec,

		    const bnd::modification* pAtpBoundMod,
		    const bnd::modification* pAdpBoundMod,
		    const bnd::modification* pPhosphorylatedMod,
		    const bnd::modification* pNoneMod,

		    twentyElevenExtrapolator* pExtrapolator) :
      generator(rMzrUnit,

		pSte20ModMol,
		ste20AtpBindingSiteNdx,
		ste20MolSpec,

		pSte11ModMol,
		rActivityMask,
		ste11MolSpec,

		pAtpBoundMod,
		pAdpBoundMod,
		pPhosphorylatedMod,
		pNoneMod,

		this,
		pExtrapolator)
    {}

    plx::omniFeature::rxnGen*
    getRxnGen(void)
    {
      return &generator;
    }
  };
}

#endif // TWENTYELEVENFAM_H
