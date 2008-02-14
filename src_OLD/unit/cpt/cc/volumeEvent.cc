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

#include "fnd/sensitivityList.hh"
#include "cpt/compartmentGraph.hh"
#include "cpt/cptReaction.hh"
#include "cpt/respondReaction.hh"
#include "cpt/cptUnit.hh"
#include "cpt/unitsMgr.hh"
#include "cpt/volumeEvent.hh"
#include "cpt/cptApp.hh"

namespace cpt
{
  volumeEvent::
  volumeEvent(const std::map<int, double>& rCompartmentToVol) :
    compartmentToVol(rCompartmentToVol)
  {}

  class updateCompartmentVolume :
    public std::unary_function<std::map<int, double>::value_type, void>
  {
    fnd::sensitivityList<cptReaction>& rAffected;
    compartmentGraph& rGraph;

  public:
    updateCompartmentVolume(fnd::sensitivityList<cptReaction>& rAffectedReactions,
			    compartmentGraph& rCompartmentGraph) :
      rAffected(rAffectedReactions),
      rGraph(rCompartmentGraph)
    {}

    void
    operator()(const argument_type& rNdxVolPair) const
    {
      int compartmentNdx = rNdxVolPair.first;
      double newVolume = rNdxVolPair.second;

      compartment* pCompartment
	= rGraph.compartments[compartmentNdx];

      pCompartment->updateVolume(newVolume,
				 rAffected);
    }
  };

  fnd::eventResult
  volumeEvent::
  happen(cptApp& rCptApp)
    throw(std::exception)
  {
    cptUnit& rCptUnit = rCptApp.getUnits()->getCptUnit();
    
    // Reset compartment volumes.
    fnd::sensitivityList<cptReaction> affectedReactions;
    std::for_each(compartmentToVol.begin(),
		  compartmentToVol.end(),
		  updateCompartmentVolume(affectedReactions,
					  rCptUnit.getCompartmentGraph()));
  
    // Reschedule affected rections.
    for_each(affectedReactions.begin(),
	     affectedReactions.end(),
	     respondReaction(rCptApp.getPropensities()));

    return fnd::go;
  }
}
