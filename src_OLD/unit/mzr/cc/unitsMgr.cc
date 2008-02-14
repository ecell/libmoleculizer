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

#include "mzr/unitsMgr.hh"

#include "mzr/mzrUnit.hh"
#include "mol/molUnit.hh"
#include "dimer/dimerUnit.hh"
#include "plex/plexUnit.hh"
#include "stoch/stochUnit.hh"
#include "ftr/ftrUnit.hh"

namespace mzr
{
  unitsMgr::
  unitsMgr(moleculizer& rMoleculizer) :
    pMzrUnit(new mzr::mzrUnit(rMoleculizer)),
    pMolUnit(new bnd::molUnit(rMoleculizer)),
    pPlexUnit(new plx::plexUnit(rMoleculizer,
				*pMzrUnit,
				*pMolUnit)),
    pStochUnit(new stoch::stochUnit(rMoleculizer,
				    *pMzrUnit)),
    pDimerUnit(new dimer::dimerUnit(rMoleculizer,
    				    *pMzrUnit,
    				    *pMolUnit,
    				    *pPlexUnit)),
    pFtrUnit(new ftr::ftrUnit(rMoleculizer,
     			      *pMzrUnit,
			      *pMolUnit,
			      *pPlexUnit))
  {
    // Note that these need to be in output order, which is slightly
    // different from linkage order: mzrUnit "plays cleanup" by parsing
    // after other units are done.

    // dimerUnit has to come before plexUnit, since the plex parser checks
    // for impossible bindings.
    //
    // ftrUnit has to come after plexUnit, since it uses omniPlexes in its
    // reaction generators.
    addEntry(pStochUnit);
    addEntry(pMolUnit);
    addEntry(pDimerUnit);
    addEntry(pPlexUnit);
    addEntry(pFtrUnit);
    addEntry(pMzrUnit);
  }

  // Class for constructing overall input capabilities of moleculizer
  // from input capabilities of each unit.
  class addCapToUnion :
    public std::unary_function<const unit*, void>
  {
    inputCapabilities& rUnionCaps;
  public:
    addCapToUnion(inputCapabilities& rUnionCapabilities) :
      rUnionCaps(rUnionCapabilities)
    {}
    
    void
    operator()(const unit* pUnit) const
    {
      rUnionCaps.addCap(pUnit->inputCap);
    }
  };
      
  void
  unitsMgr::
  unionInputCaps(inputCapabilities& rUnion)
  {
    for_each(begin(),
	     end(),
	     addCapToUnion(rUnion));
  }
}
