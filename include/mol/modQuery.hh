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

#ifndef MODQUERY_H
#define MODQUERY_H

#include "plex/plexQuery.hh"
#include "mol/molState.hh"
#include "mol/modMixin.hh"

namespace bnd
{
  class modMixinStateQuery : public plx::plexQuery
  {
    plx::plexMolSpec molSpec;

    // The index of the modification site to look at.
    // This is generated from the name when the query is created.
    int modNdx;
    // Pointer to the modification to see at the modification site.
    // This is generated (by getMod) when the query is created.
    const modification* pMod;

  public:

    modMixinStateQuery(const plx::plexMolSpec& rPlexMolSpec,
		       int modificationIndex,
		       const modification* pModToSee) :
      molSpec(rPlexMolSpec),
      modNdx(modificationIndex),
      pMod(pModToSee)
    {}

    bool
    operator()(const plx::plexParam& rParam) const
    {
      // Get the modMixin part of the state of the mol.
      const modStateMixin* pModMixin
	= dynamic_cast<const modStateMixin*>(rParam.molParams[molSpec]);

      // Complain if the mol state does not have a modMixin part.
      if(! pModMixin) throw modMolQueryTargetXcpt(pMod,
						  modNdx);
      return (*pModMixin)[modNdx] == pMod;
    }

    // Why did I not overload operator() ?
    bool
    applyTracked(const plx::plexParam& rPlexParam,
		 const plx::subPlexSpec& rSubPlexSpec) const
    {
      // Map the mol index through the injection of the omniplex into
      // the complex.
      int molNdx = rSubPlexSpec.getInjection().forward.molMap[molSpec];

      // Get the state of the mol at the mapped mol index.
      const molState* pMolState = rPlexParam.molParams[molNdx];

      // Get the modMixin part of the state of the mol.
      const modStateMixin* pModMixin
	= dynamic_cast<const modStateMixin*>(pMolState);

      if(! pModMixin) throw modMolQueryTargetXcpt(pMod,
						  modNdx);
      return (*pModMixin)[modNdx] == pMod;
    }
  };
}

#endif // MODQUERY_H
