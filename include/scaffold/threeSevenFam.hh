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

#ifndef THREESEVENFAM_H
#define THREESEVENFAM_H

#include "mzr/reactionFamily.hh"
#include "scaffold/threeSevenRxnGen.hh"

namespace scaf
{
  class threeSevenFam : public mzr::reactionFamily
  {
    threeSevenRxnGen generator;

  public:

    threeSevenFam(mzr::mzrUnit& rMzrUnit,

		  bnd::modMol* pFus3ModMol,
		  int fus3AtpModSiteNdx,
		  plx::plexMolSpec fus3MolSpec,
		  int minFus3ActivePhosCount,

		  bnd::modMol* pSte7ModMol,
		  plx::plexMolSpec ste7MolSpec,
		  std::vector<bool>& rActivityMask,

		  const bnd::modification* pAtpBoundMod,
		  const bnd::modification* pAdpBoundMod,
		  const bnd::modification* pPhosphorylatedMod,
		  const bnd::modification* pNoneMod,

		  threeSevenExtrapolator* pExtrapolator) :
      generator(rMzrUnit,

		pFus3ModMol,
		fus3AtpModSiteNdx,
		fus3MolSpec,
		minFus3ActivePhosCount,

		pSte7ModMol,
		ste7MolSpec,
		rActivityMask,

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

#endif // THREESEVENFAM_H
