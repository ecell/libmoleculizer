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

#ifndef NUCEXFAM_H
#define NUCEXFAM_H

#include "mzr/reactionFamily.hh"
#include "nucEx/nucExRxnGen.hh"

namespace nucEx
{
  // Family of nucleotide-exchange reactions which can be used to handle
  // Gpa1's nucleotide-exchange reaction.

  class nucExFam : public mzr::reactionFamily
  {
    nucExRxnGen generator;
  public:
    nucExFam(bnd::modMol* pTargetMol,
	     int nucModSiteNdx,
	     const bnd::modification* pNoneModification,
	     const bnd::modification* pBoundModification,
	     plx::plexFamily* pEnablingSubcomplex,
	     plx::plexMolSpec enablingMolSpec,
	     stoch::stochSpecies* pNucleotideSpecies,
	     mzr::mzrUnit& rMzrUnit,
	     nucExBindExtrapolator* pEnabledBindExtrapolator,
	     nucExUnbindExtrapolator* pEnabledUnbindExtrapolator,
	     nucExBindExtrapolator* pPlainBindExtrapolator,
	     nucExUnbindExtrapolator* pPlainUnbindExtrapolator) :
      generator(pTargetMol,
		nucModSiteNdx,
		pNoneModification,
		pBoundModification,
		pEnablingSubcomplex,
		enablingMolSpec,
		pNucleotideSpecies,
		this,
		rMzrUnit,
		pEnabledBindExtrapolator,
		pEnabledUnbindExtrapolator,
		pPlainBindExtrapolator,
		pPlainUnbindExtrapolator)
    {}

    plx::molInContext::rxnGen*
    getRxnGen(void)
    {
      return &generator;
    }
  };
}

#endif // NUCEXFAM_H
