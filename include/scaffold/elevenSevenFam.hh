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

#ifndef ELEVENSEVENFAM_H
#define ELEVENSEVENFAM_H

#include "mzr/reactionFamily.hh"
#include "scaffold/elevenSevenRxnGen.hh"

namespace scaf
{
  class elevenSevenFam : public mzr::reactionFamily
  {
    elevenSevenRxnGen generator;

  public:
    elevenSevenFam(mzr::mzrUnit& rMzrUnit,

		   bnd::modMol* pSte11ModMol,
		   int ste11AtpModSiteNdx,
		   plx::plexMolSpec ste11MolSpec,
		   int elevenActivePhosCount,

		   bnd::modMol* pSte7ModMol,
		   plx::plexMolSpec ste7MolSpec,
		   int ste7TargetModSiteNdx,

		   const bnd::modification* pAtpBoundMod,
		   const bnd::modification* pAdpBoundMod,
		   const bnd::modification* pPhosphorylatedMod,
		   const bnd::modification* pNoneMod,

		   elevenSevenExtrapolator* pExtrapolator) :
      generator(rMzrUnit,

		pSte11ModMol,
		ste11AtpModSiteNdx,
		ste11MolSpec,
		elevenActivePhosCount,

		pSte7ModMol,
		ste7MolSpec,
		ste7TargetModSiteNdx,

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

#endif // ELEVENSEVENFAM_H
