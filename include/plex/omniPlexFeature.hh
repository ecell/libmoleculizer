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

#ifndef OMNIPLEXFEATURE_H
#define OMNIPLEXFEATURE_H

/*! \file omniPlexFeature.hh
  \ingroup omniGroup
  \brief Defines subcomplex feature. */

#include "mzr/feature.hh"
#include "mzr/dumpable.hh"
#include "plex/plexSpecies.hh"
#include "plex/subPlexSpec.hh"

namespace plx
{
  class omniPlexFeature :
    public mzr::feature<plx::plexSpecies, plx::subPlexSpec>
  {
    // If a dumpable is really attached to this feature, then
    // thiw points to it; otherwise null.
    mzr::multiSpeciesDumpable* pDumpable;
    
  public:
    omniPlexFeature(void) :
      pDumpable(0)
    {}

    // Generate reactions and add new species to the dumpable,
    // if any.
    void
    notifyNew(plx::plexSpecies* pNewSpecies,
	      const plx::subPlexSpec& rSpec);

    // To "turn on" dumping of the species in this omniplex.
    void
    setDumpable(mzr::multiSpeciesDumpable* ptrDumpable)
    {
      pDumpable = ptrDumpable;
    }
  };
}

#endif
