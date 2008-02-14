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

#include "cpt/cptReaction.hh"
#include "cpt/respondReaction.hh"
#include "cpt/growEvent.hh"
#include "cpt/unitsMgr.hh"
#include "cpt/cptUnit.hh"
#include "cpt/cptApp.hh"

namespace cpt
{
  growEvent::
  growEvent(double growthFactor,
	    double schedulingPeriod) :
    factor(growthFactor),
    period(schedulingPeriod)
  {}

  class growCompartment :
    public std::unary_function<compartment*, void>
  {
    double factor;
    fnd::sensitivityList<cptReaction>& rAffected;

  public:
    growCompartment(double growthFactor,
		    fnd::sensitivityList<cptReaction>& rAffectedReactions) :
      factor(growthFactor),
      rAffected(rAffectedReactions)
    {
    }

    void
    operator()(compartment* pCompartment) const
    {
      double currentVolume = pCompartment->getVolume();
      
      pCompartment->updateVolume(factor * currentVolume,
				 rAffected);
    }
  };

  fnd::eventResult
  growEvent::
  happen(cptApp& rCptApp)
    throw(std::exception)
  {
    const compartmentGraph& rGraph
      = rCptApp.getUnits()->getCptUnit().getCompartmentGraph();

    // Grow all the compartments by the same factor.
    fnd::sensitivityList<cptReaction> affectedReactions;
    std::for_each(rGraph.compartments.begin(),
		  rGraph.compartments.end(),
		  growCompartment(factor,
				  affectedReactions));

    // Update the propensities of reactions.
    for_each(affectedReactions.begin(),
	     affectedReactions.end(),
	     respondReaction(rCptApp.getPropensities()));

    // Reschedule this event for later.
    rCptApp.getQueue().scheduleEvent(this,
				     rCptApp.getSimTime() + period);

    return fnd::go;
  }
}
