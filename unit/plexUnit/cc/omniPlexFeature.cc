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

#include "plex/plexQuery.hh"
#include "plex/omniPlex.hh"
#include "plex/omniPlexFeature.hh"

namespace plx
{
  void
  omniPlexFeature::notifyNew(plx::plexSpecies* pNewSpecies,
			     const plx::subPlexSpec& rSpec)
  {
    // Get the omniplex's state query.
    //
    // This seems the long way around, since this feature will be attached
    // directly to the omniplex we get here by rSpec.getOmni().
    const andPlexQueries* pQuery = rSpec.getOmni()->getStateQuery();

    // Create new reactions and add species to dumpable if species state
    // satisfies the state query.  We apply the state query through the
    // imbedding of the omniplex into the new species.
    if(pQuery->applyTracked(pNewSpecies->getParam(),
			    rSpec))
      {
	// Create the new reactions.
	mzr::feature<plx::plexSpecies, plx::subPlexSpec>::
	  notifyNew(pNewSpecies,
		    rSpec);

	// If there is a dumpable, notify it of the new species.
	if(pDumpable) pDumpable->addSpecies(pNewSpecies);
      }
  }
}
