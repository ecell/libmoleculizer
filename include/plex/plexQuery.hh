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

#ifndef PLEXQUERY_H
#define PLEXQUERY_H

#include <algorithm>
#include "mzr/species.hh"
#include "mzr/util.hh"
#include "mzr/paramDumpable.hh"
#include "plex/subPlexSpec.hh"
#include "plex/prm.hh"

namespace plx
{
  class plexUnit;

  /*! \file plexQuery.hh
    \ingroup plexSpeciesGroup
    \brief Defines logical combinators for plex dump queries. */

  /*! \ingroup plexSpeciesGroup
    \brief Base class for queries about complexes.

    Adds to plexSpecies::paramQuery the need to map plexParam argument
    through a plex isomorphism. This way of applying the query is used
    in omniplex dumps.
  */
  class plexQuery : public mzr::paramQuery<plexParam>
  {
  public:
    plexQuery(plexUnit& rPlexUnit);

    // Seem to be having some trouble overloading operator() to handle
    // this case.
    virtual bool applyTracked(const plexParam& rPlexParam,
			      const subPlexSpec& rSpec) const = 0;
  };

  /*! \ingroup plexSpeciesGroup
    \brief AND logical combinator for plex dump queries.

    The command only supports AND of a list of canned atomic queries. */
  class andPlexQueries : public plexQuery
  {
    // Pointers to the queries that this query AND's together.
    // This query is NOT responsible for memory-managing these
    // "subordinate" queries; since plexQueries can refer
    // to one another, they are all registered in the plexUnit.
    std::vector<plexQuery*> queries;

  public:
    andPlexQueries(plexUnit& rPlexUnit) :
      plexQuery(rPlexUnit)
    {}

    void addQuery(plexQuery* pQuery)
    {
      queries.push_back(pQuery);
    }

    bool operator()(const plexParam& rParam) const;

    bool applyTracked(const plexParam& rParam,
		      const subPlexSpec& rSpec) const;
  };
}

#endif
