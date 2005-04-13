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

#include "mol/molQuery.hh"
#include "mol/molUnit.hh"

namespace bnd
{
  molQuery::molQuery(molUnit& rMolUnit)
  {
    rMolUnit.addMolQuery(this);
  }
  
  // This is accessible via adaptors, I think.
  class applyNotMolQuery :
    public std::unary_function<const molQuery*, bool>
  {
    const molParam& rParam;
  public:
    applyNotMolQuery(const molParam& rMolParam) :
      rParam(rMolParam)
    {}

    bool
    operator()(const molQuery* pQuery) const
    {
      const molQuery& rQuery = *pQuery;
      return ! rQuery(rParam);
    }
  };

  bool
  andMolQueries::
  operator()(const molParam& rParam) const
  {
    return (queries.end() == find_if(queries.begin(),
				     queries.end(),
				     applyNotMolQuery(rParam)));
  }
}
