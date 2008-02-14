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

#ifndef CDM_DECOMPRXNGEN_H
#define CDM_DECOMPRXNGEN_H

#include "clx/clxUnit.hh"
#include "cdm/decomposeExtrap.hh"

namespace cdm
{
  /*! \ingroup decompGroup
    \brief Reaction generator decompositions. */
  class decompRxnGen :
    public fnd::rxnGen<cpx::cxBinding<clx::cptPlexSpecies, clx::cptPlexFamily> >
  {
    utl::autoVector<cpt::globalReaction>* pFamily;
    cpt::cptApp& rCptApp;
    cpt::cptUnit& rCptUnit;
    clx::clxUnit& rClxUnit;
    decomposeExtrapolator* pExtrap;
  
  public:
    // Note that this reaction generator memory manages the rate extrapolator.
    decompRxnGen(utl::autoVector<cpt::globalReaction>* pDecompFamily,
		 cpt::cptApp& refCptApp,
		 cpt::cptUnit& refCptUnit,
		 clx::clxUnit& refClxUnit,
		 decomposeExtrapolator* pExtrapolator) :
      pFamily(pDecompFamily),
      rCptApp(refCptApp),
      rCptUnit(refCptUnit),
      rClxUnit(refClxUnit),
      pExtrap(pExtrapolator)
    {}
    
    ~decompRxnGen(void)
    {
      delete pExtrap;
    }

    void
    respond(const fnd::featureStimulus<cpx::cxBinding<clx::cptPlexSpecies, clx::cptPlexFamily> >& rStimulus);

  };
}

#endif // CDM_DECOMPRXNGEN_H
