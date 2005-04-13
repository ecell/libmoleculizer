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

#include "mol/modQuery.hh"
#include "mol/modMixin.hh"

namespace bnd
{
  bool
  modMixinMolQuery::
  operator()(const molParam& rParam) const
  {
    // Get the modMixin part of the state of the mol.
    const modStateMixin* pModMixin
      = dynamic_cast<const modStateMixin*>(rParam);

    // Complain if the mol state does not have a modMixin part.
    if(! pModMixin) throw modMolQueryTargetXcpt(pMod,
						modNdx);
    return (*pModMixin)[modNdx] == pMod;
  }
  
  bool
  modMixinStateQuery::
  operator()(const plx::plexParam& rParam) const
  {
    const modMixinMolQuery& rMolQuery = *pMolQuery;
    return rMolQuery(rParam.molParams[molSpec]);
  }

  bool
  modMixinStateQuery::
  applyTracked(const plx::plexParam& rPlexParam,
	       const plx::subPlexSpec& rSubPlexSpec) const
  {
    // Map the mol index through the injection of the omniplex into
    // the complex.
    int molNdx = rSubPlexSpec.getInjection().forward.molMap[molSpec];

    // Test the state of the mol at the mapped mol index.
    const modMixinMolQuery& rMolQuery = *pMolQuery;
    return rMolQuery(rPlexParam.molParams[molNdx]);
  }
}
