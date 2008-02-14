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

#ifndef CFT_UNIMOLFAM_H
#define CFT_UNIMOLFAM_H

#include "utl/autoVector.hh"
#include "cpt/globalSpecies.hh"
#include "cpt/globalReaction.hh"
#include "cpt/cptUnit.hh"
#include "clx/clxUnit.hh"
#include "cml/cptModMol.hh"
#include "cft/uniMolGen.hh"

namespace cft
{
  class uniMolFam :
    public utl::autoVector<cpt::globalReaction>
  {
    uniMolRxnGen rxnGen;

  public:
    uniMolFam(cpt::cptApp& refCptApp,
	      cpt::cptUnit& refCptUnit,
	      clx::clxUnit& refClxUnit,
	      cml::cptModMol* pEnablingModMol,
	      cpx::andMolStateQueries* pAndMolQueries,
	      const std::vector<molModExchange>& rModExchanges,
	      cpt::globalSpecies* pAuxReactant,
	      cpt::globalSpecies* pAuxProduct,
	      const uniMolExtrapolator* pUniMolExtrapolator) :
      rxnGen(refCptApp,
	     refCptUnit,
	     refClxUnit,
	     pEnablingModMol,
	     pAndMolQueries,
	     rModExchanges,
	     pAuxReactant,
	     pAuxProduct,
	     this,
	     pUniMolExtrapolator)
    {}

    uniMolRxnGen*
    getRxnGen(void)
    {
      return &rxnGen;
    }
  };
}

#endif // CFT_UNIMOLFAM_H
