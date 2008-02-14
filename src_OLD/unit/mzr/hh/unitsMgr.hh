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

#ifndef UNITSMGR_H
#define UNITSMGR_H

#include "utl/autoVector.hh"
#include "mzr/unit.hh"

namespace mzr
{
  class mzrUnit;
}
namespace bnd
{
  class molUnit;
}
namespace plx
{
  class plexUnit;
}
namespace dimer
{
  class dimerUnit;
}
namespace stoch
{
  class stochUnit;
}
namespace ftr
{
  class ftrUnit;
}

namespace mzr
{
  class moleculizer;
  
  // As a vector, the units should be ordered in output order, so that if each
  // unit arranges its elements correctly (under various heads in the output
  // document) and the units are ordered thus, then all the elements in the
  // output doc should appear in the order demanded by the schema.
  class unitsMgr :
    public utl::autoVector<unit>
  {
  public:

    mzr::mzrUnit* pMzrUnit;
    bnd::molUnit* pMolUnit;
    plx::plexUnit* pPlexUnit;
    stoch::stochUnit* pStochUnit;
    dimer::dimerUnit* pDimerUnit;
    ftr::ftrUnit* pFtrUnit;

    unitsMgr(moleculizer& rMoleculizer);

    // Prepares the overall input capabilities of moleculizer from
    // the several input capabilities of the units.
    void
    unionInputCaps(inputCapabilities& rUnion);
  };
}

#endif // UNITSMGR_H
