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

#ifndef CPT_CREATEEVENT_H
#define CPT_CREATEEVENT_H

#include "cpt/cptEvent.hh"
#include "cpt/globalSpecies.hh"
#include "cpt/compartmentGraph.hh"

namespace cpt
{
  // Event that creates molecules (at the beginning or) during
  // a simulation.
  class createEvent :
    public cptEvent
  {
    globalSpecies* pSpeciesToCreate;

    // This bodes ill for the idea of reorganizing compartments on the fly.  A
    // map from compartment pointers to populations might work better,
    // particularly if one only subdivides compartments and never merges
    // them. With that restriction, each compartment as an entity could
    // continue to exist (reasonably sensibly) throughout the simulation.

    // This vector's indexing should match that of the compartment
    // vector in the compartment graph.
    std::vector<int> compartmentPops;

  public:
    // In order to create specified numbers of molecules in the different
    // compartments.
    //
    // The pops vector should be indexed like the compartments vector
    // in the compartment graph.
    createEvent(globalSpecies* pSpecies,
		const std::vector<int>& rCompartmentPops);

    // In order to create the same number of molecules in every compartment.
    createEvent(globalSpecies* pSpecies,
		int defaultCompartmentPop);

    fnd::eventResult
    happen(cptApp& rApp)
      throw(std::exception);
  };
}

#endif // CPT_CREATEEVENT_H
