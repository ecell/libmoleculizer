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
#include "plex/plexUnit.hh"

namespace plx
{
  plexQuery::plexQuery(plexUnit& rPlexUnit)
  {
    rPlexUnit.addPlexQuery(this);
  }

  /*! \brief Auxiliary function class to help do AND of a number of
    plexQueries. */
  class applyNotQuery : public std::unary_function<const plexQuery*, bool>
  {
    const plexParam& rParam;
  public:
    applyNotQuery(const plexParam& rPlexParam) :
      rParam(rPlexParam)
    {}

    bool
    operator()(const plexQuery* pQuery) const
    {
      const plexQuery& rQuery = *pQuery;
      return ! rQuery(rParam);
    }
  };

  bool
  andPlexQueries::
  operator()(const plexParam& rParam) const
  {
    // The AND of the queries is true if we can get all the way to
    // the end of the vector of queries without finding one that is not
    // satisfied by the parameter.
    return (queries.end() ==  find_if(queries.begin(),
				      queries.end(),
				      applyNotQuery(rParam)));
  }

  /*! \brief Auxiliary function class to help do AND of a number of
    plexQueries in the case that an isomorphism must be applied to
    the plexParam argument. */
  class applyNotTrackedQuery :
    public std::unary_function<const plexQuery*, bool>
  {
    const plexParam& rParam;
    const subPlexSpec& rSpec;
  public:
    applyNotTrackedQuery(const plexParam& rPlexParam,
			 const subPlexSpec& rSubPlexSpec) :
      rParam(rPlexParam),
      rSpec(rSubPlexSpec)
    {}

    bool
    operator()(const plexQuery* pQuery) const
    {
      const plexQuery& rQuery = *pQuery;
      return ! rQuery.applyTracked(rParam,
				   rSpec);
    }
  };

  bool
  andPlexQueries::
  applyTracked(const plexParam& rParam,
	       const subPlexSpec& rSpec) const
  {
    return (queries.end() == find_if(queries.begin(),
				     queries.end(),
				     applyNotTrackedQuery(rParam,
							  rSpec)));
  }
}
