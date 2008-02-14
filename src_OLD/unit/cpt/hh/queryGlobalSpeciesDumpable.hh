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

#ifndef CPT_QUERYGLOBALSPECIESDUMPABLE_H
#define CPT_QUERYGLOBALSPECIESDUMPABLE_H

#include "fnd/query.hh"
#include "cpt/multiGlobalSpeciesDumpable.hh"
#include "cpt/compartmentGraph.hh"

namespace cpt
{
  // These dumpable classes apparently have to be reimplemented for use with
  // globalSpecies, even though the member functions that don't work with
  // globalSpecies don't really need to be instantiated.  g++ basically
  // instantiates these template member functions as part of its way of
  // handling the overall linkage problem with templates.  (For example, the
  // diagnostic does not say where the code was instantiated, as it always
  // does when the template MUST be instantiated.

  // This one has to be a template, because the query that controls
  // which species get dumped depends on the precise type of species.
  template<class globalSpeciesT>
  class queryGlobalSpeciesDumpable :
    public multiGlobalSpeciesDumpable<globalSpeciesT>
  {
  public:
    typedef globalSpeciesT globalSpeciesType;
    
  protected:
    fnd::query<globalSpeciesType>& rQuery;
    
  public:
    queryGlobalSpeciesDumpable(const std::string& rName,
			       const compartmentGraph& rCompartmentGraph,
			       fnd::query<globalSpeciesType>& rGlobalSpeciesQuery) :
      multiGlobalSpeciesDumpable<globalSpeciesType>(rName,
						    rCompartmentGraph),
      rQuery(rGlobalSpeciesQuery)
    {}

    void
    respond(const fnd::newSpeciesStimulus<globalSpeciesType>& rStimulus)
    {
      const globalSpeciesType* pNewSpecies
	= rStimulus.getSpecies();
      
      if(rQuery(*pNewSpecies))
	{
	  this->dumpedSpecies.push_back(pNewSpecies);
	}
    }
  };
}

#endif // CPT_QUERYGLOBALSPECIESDUMPABLE_H
