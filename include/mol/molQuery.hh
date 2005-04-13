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

#ifndef MOLQUERY_H
#define MOLQUERY_H

#include "mol/molState.hh"
#include "mzr/paramDumpable.hh"

namespace bnd
{
  class molUnit;

  class molQuery :
    public mzr::paramQuery<molParam>
  {
  public:
    molQuery(molUnit& rMolUnit);
  };

  class andMolQueries : public molQuery
  {
    // These have to be stored by reference, since molQuery is just a base
    // class.  I intend for molQueries to be managed by the mol unit.
    std::vector<const molQuery*> queries;
    
  public:
    andMolQueries(molUnit& rMolUnit) :
      molQuery(rMolUnit)
    {}
    
    // Note that this andMolQueries does not memory-manage queries added to
    // it.
    void
    addQuery(const molQuery* pQuery)
    {
      queries.push_back(pQuery);
    }

    bool
    operator()(const molParam& rParam) const;
  };
}

#endif // MOLQUERY_H
