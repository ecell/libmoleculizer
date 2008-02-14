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

#ifndef CLX_CPTOMNISTRUCTUREQUERIES_H
#define CLX_CPTOMNISTRUCTUREQUERIES_H

#include "fnd/query.hh"
#include "cpx/omniStructureQuery.hh"
#include "clx/cptPlex.hh"

namespace clx
{
  // Query on the structure of a plex that putatively matches
  // an omniPlex.  Only one type of structural query exists now:
  // is a site that is free on the omniplex actually free in the
  // complex where the omniplex is found. (omniFreeSiteQuery below).
  typedef
  fnd::query<cpx::omniStructureQueryArg<cptPlex> >
  cptOmniStructureQuery;

  // The only kind of structural query used by omniplexes at this time.
  typedef
  cpx::omniFreeSiteQuery<cptPlex>
  cptOmniFreeSiteQuery;

  // Conjunction of structural queries.  Each omniPlex has one of these
  // to "gate" putatively matching plexSpecies.
  typedef
  fnd::andQueries<cptOmniStructureQuery>
  cptOmniStructureQueries;
}

#endif // CLX_CPTOMNISTRUCTUREQUERIES_H
