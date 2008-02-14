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

#ifndef CPTAPP_H
#define CPTAPP_H

#include "utl/gsl.hh"
#include "utl/domJob.hh"
#include "utl/autoCatalog.hh"
#include "utl/autoVector.hh"
#include "cpt/globalSpecies.hh"
#include "cpt/globalReaction.hh"
#include "cpt/propensityDistro.hh"
#include "cpt/eventQueue.hh"
#include "cpt/inputCapabilities.hh"

namespace cpt
{
  class unitsMgr;

  class cptUnit;

  class cptApp :
    public utl::dom::domBatchApp
  {
    // Assembles the input capabilities from the units.
    // 
    // This was all intended to help make sure that the constructs present in
    // a given input file could actually be handled by the units that happened
    // to be in this program.
    void
    constructorPrelude(void);

    // Verifies input elements can be handled by units.  Invokes each unit's
    // parsing routine.
    void
    constructorCore(xmlpp::Element* pRootElement,
		    xmlpp::Element* pModelElement,
		    xmlpp::Element* pStreamsElement,
		    xmlpp::Element* pEventsElement)
      throw(std::exception);

    // Things that are legal in cpt input.  This is assembled from
    // information provided by the linked-in units.
    inputCapabilities inputCap;
    
    // The libraries out of which cpt is built.
    unitsMgr* pUnits;

    // For easier/quicker access to the mandatory unit.
    cptUnit& rCptUnit;

    // Run parameters that come from the input file.
    double cycleTime;
    double reactionsPerCycle;
    double epsilon;

    // The current simulation time.
    double now;

    propensityDistro propensities;

    // This is for user-scheduled events, such as dumping output, dumping
    // state, or stopping.  Reactions in cpt are not events in this sense.
    eventQueue theQueue;

    // Inner core of "run." This routine exists to enable easy exit
    // from the simulation cycle (by returning from this routine.)
    void
    simulationCycle(utl::gsl::multinomialSampler& rSampler,
		    int& rCycleCount,
		    int& rEventCount);

  public:
    cptApp(int argCount,
	   char** argVector,
	   xmlpp::Document* pDoc)
      throw(std::exception);

    // We need to delete the units when we go away.
    ~cptApp(void);

    void
    parseCommandLine(int argc,
		     char* argv[]);

    // For access to all the units.
    unitsMgr*
    getUnits(void)
    {
      return pUnits;
    }

    // Quicker access to the only mandatory unit.
    cptUnit&
    getCptUnit(void)
    {
      return rCptUnit;
    }

    // Gets the maximum time step, parsed out of the input file.
    double
    getCycleTime(void)
    {
      return cycleTime;
    }

    // Get "target" number of reaction events per cycle.
    // This is really used, together with the total propensity of all the
    // reactions, to limit the length of the cycle.
    double
    getReactionsPerCycle(void)
    {
      return reactionsPerCycle;
    }

    double
    getEpsilon(void)
    {
      return epsilon;
    }

    propensityDistro&
    getPropensities(void)
    {
      return propensities;
    }

    // Access to the user event queue, which holds events scheduled by the
    // user, including stop event, dump output events, dump state event, etc.,
    // but not reactions, as in Moleculizer.
    eventQueue&
    getQueue(void)
    {
      return theQueue;
    }

    const eventQueue&
    getQueue(void) const
    {
      return theQueue;
    }

    double
    getSimTime(void) const
    {
      return now;
    }

    void
    setSimTime(double newSimTime)
    {
      now = newSimTime;
    }

    // Generate a state dump.
    xmlpp::Document*
    makeDomOutput(void)
      throw(std::exception);

    // Run the simulation.
    int
    run(void)
      throw(std::exception);
  };

}

#endif // CPTAPP_H
