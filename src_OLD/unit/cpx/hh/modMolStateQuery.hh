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

#ifndef CPX_MODMOLQUERY_H
#define CPX_MODMOLQUERY_H

#include "cpx/molStateQuery.hh"
#include "cpx/modStateMixin.hh"
#include "cpx/badModMolStateXcpt.hh"

namespace cpx
{
  // Query to test if a modMol has a particular modification at a particular
  // modification site.
  class modMolStateQuery :
    public molStateQuery
  {
    int modNdx;
    const modification* pMod;

  public:
    modMolStateQuery(int modificationIndex,
		     const modification* pModToSee) :
      modNdx(modificationIndex),
      pMod(pModToSee)
    {}

    bool
    operator()(const molParam& rParam) const
    {
      const modStateMixin* pModStateMixin
	= dynamic_cast<const modStateMixin*>(rParam);

      if(! pModStateMixin) throw badModMolStateXcpt();

      return (*pModStateMixin)[modNdx] == pMod;
    }
  };
}

#endif // CPX_MODMOLQUERY_H
