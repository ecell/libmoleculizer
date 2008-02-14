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

#ifndef CPT_UNITSMGR_H
#define CPT_UNITSMGR_H

#include "utl/autoVector.hh"
#include "cpt/unit.hh"

// Compartmentl analogue of mzrUnit.
namespace cpt
{
  class cptUnit;
}

// Compartmental version of stochUnit
namespace cst
{
  class cstUnit;
}

// Compartmental version of molUnit.
namespace cml
{
  class cmlUnit;
}

// Compartmental version of plexUnit.
namespace clx
{
  class clxUnit;
}

// Compartmental version of dimerUnit.
namespace cdm
{
  class cdmUnit;
}

// Compartmental version of ftrUnit.
namespace cft
{
  class cftUnit;
}

namespace cpt
{
  class cptApp;
  
  // As a vector, the units should be ordered in output order, so that if each
  // unit arranges its elements correctly (under various heads in the output
  // document) and the units are ordered thus, then all the elements in the
  // output doc should appear in the order demanded by the schema.
  class unitsMgr :
    public utl::autoVector<unit>
  {
    cpt::cptUnit* pCptUnit;
    cst::cstUnit* pCstUnit;
    cml::cmlUnit* pCmlUnit;
    clx::clxUnit* pClxUnit;
    cdm::cdmUnit* pCdmUnit;
    cft::cftUnit* pCftUnit;

  public:

    unitsMgr(cptApp& rCptApp);

    cpt::cptUnit& 
    getCptUnit(void)
    {
      return *pCptUnit;
    }

    cst::cstUnit&
    getCstUnit(void)
    {
      return *pCstUnit;
    }

    cml::cmlUnit&
    getCmlUnit(void)
    {
      return *pCmlUnit;
    }

    clx::clxUnit&
    getClxUnit(void)
    {
      return *pClxUnit;
    }

    cft::cftUnit&
    getCftUnit(void)
    {
      return *pCftUnit;
    }

    // Prepares the overall input capabilities of moleculizer from
    // the several input capabilities of the units.
    void
    unionInputCaps(inputCapabilities& rUnion);
  };
}

#endif // CPT_UNITSMGR_H
