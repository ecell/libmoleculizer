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

#include "mzr/mzrUnit.hh"
#include "mzr/unitsMgr.hh"
#include "mzr/moleculizer.hh"
#include "mzr/tabDumpEvent.hh"
#include "mzr/dupSpeciesNameXcpt.hh"
#include "mzr/unkSpeciesXcpt.hh"
#include "mzr/dupDumpableNameXcpt.hh"
#include "mzr/unkDumpableXcpt.hh"
#include "mzr/clockDumpable.hh"
#include "mzr/secondsDumpable.hh"
#include "mzr/reactEventCountDumpable.hh"
#include "mzr/reactCountDumpable.hh"
#include "mzr/speciesCountDumpable.hh"
#include "mzr/simTimeDumpable.hh"
#include "mzr/volumeDumpable.hh"

namespace mzr
{
  void
  mzrUnit::
  mustAddSpecies(const std::string& rSpeciesName,
		 mzrSpecies* pSpecies,
		 xmlpp::Node* pRequestingNode)
    throw(utl::xcpt)
  {
    if(! addSpecies(rSpeciesName,
		    pSpecies))
      throw dupSpeciesNameXcpt(rSpeciesName,
			       pRequestingNode);
  }

  mzrSpecies*
  mzrUnit::
  mustFindSpecies(const std::string& rSpeciesName,
		  xmlpp::Node* pRequestingNode) const
    throw(utl::xcpt)
  {
    mzrSpecies* pSpecies = findSpecies(rSpeciesName);

    if(! pSpecies)
      throw unkSpeciesXcpt(rSpeciesName,
			   pRequestingNode);
    return pSpecies;
  }

  void
  mzrUnit::
  mustAddDumpable(fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable,
		  xmlpp::Node* pRequestingNode)
    throw(utl::xcpt)
  {
    if(! addDumpable(pDumpable))
      throw dupDumpableNameXcpt(pDumpable->getName(),
				pRequestingNode);
  }
  
  fnd::dumpable<fnd::basicDumpable::dumpArg>*
  mzrUnit::
  mustFindDumpable(const std::string& rDumpableName,
		   xmlpp::Node* pRequestingNode) const
    throw(utl::xcpt)
  {
    fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable
      = findDumpable(rDumpableName);

    if(! pDumpable)
      throw unkDumpableXcpt(rDumpableName,
			    pRequestingNode);

    return pDumpable;
  }

  mzrUnit::
  mzrUnit(moleculizer& rMoleculizer) :
    unit("mzr",
	 rMoleculizer),
    generateDepth(0),
    generateOption(true),
    generateOk(true),
    startSeconds(time(0)),
    secondsLimit(0)
  {
    // Initialize the vector of global state variables, with which each
    // reaction is initialized.  At this point, there is only one
    // global state variable, the volume.
    globalVars.push_back(&theMolarFactor);
    
    // Model elements whose contents are parsed by some unit
    // or another, as determined by moleculizer::parseDomInput.
    inputCap.addModelContentName(eltName::reactionGens);
    inputCap.addModelContentName(eltName::explicitSpecies);

    // Model elements that this unit actually processes.
    inputCap.addModelContentName(eltName::explicitReactions);
    inputCap.addModelContentName(eltName::volume);

    // This unit is not responsible for any reaction generators
    // or species streams.

    // Events for which this unit is responsible.
    inputCap.addEventsContentName(eltName::createEvent);
    inputCap.addEventsContentName(eltName::dumpStateEvent);
    inputCap.addEventsContentName(eltName::stopEvent);

    // Add the standard, built-in dumpables.
    addDumpable(new clockDumpable());
    addDumpable(new secondsDumpable(*this));
    addDumpable(new reactEventCountDumpable());
    addDumpable(new reactCountDumpable());
    addDumpable(new speciesCountDumpable());
    addDumpable(new simTimeDumpable(rMolzer));
    addDumpable(new volumeDumpable(*this));
  }

  // Cause tabDumpEvent to emit its header to its file and does
  // the initial scheduling of the tabDumpEvent.
  class initScheduleTabDumpEvent :
    public std::unary_function<tabDumpEvent*, void>
  {
    moleculizer& rMolzer;
  public:
    initScheduleTabDumpEvent(moleculizer& rMoleculizer) :
      rMolzer(rMoleculizer)
    {}

    void
    operator()(tabDumpEvent* pTabDumpEvent) const
    {
      pTabDumpEvent->init();
      
      rMolzer.eventQ.scheduleEvent(pTabDumpEvent,
				   rMolzer.eventQ.getSimTime());
    }
  };

  void
  mzrUnit::prepareToRun(xmlpp::Element* pRootElt,
			xmlpp::Element* pModelElt,
			xmlpp::Element* pStreamsElt,
			xmlpp::Element* pEventsElt) throw(std::exception)
  {
    // Have each tabDumpEvent emit its header line, and schedule it for time 0,
    // after which it will reschedule itself periodically.  .
    std::for_each(tabDumpEvents.begin(),
		  tabDumpEvents.end(),
		  initScheduleTabDumpEvent(rMolzer));

    // Shut off reaction generation if called for by program option.
    setGenerateOk(getGenerateOption());
}

  // Does the initial scheduling of tabDumpEvent without emitting
  // headers.  This is so that when a simulation is restarted from
  // a state dump, the output files continue instead of starting
  // at the time of the state dump.
  //
  // The tabDumpEvent is scheduled at the current simulation time,
  // which should be the time at which the state dump happened.
  class continueTabDumpEvent :
    public std::unary_function<tabDumpEvent*, void>
  {
    moleculizer& rMolzer;
  public:
    continueTabDumpEvent(moleculizer& rMoleculizer) :
      rMolzer(rMoleculizer)
    {
    }

    void
    operator()(tabDumpEvent* pTabDumpEvent) const
    {
      pTabDumpEvent->init();
      
      // Schedule the dump event at the next time that it would
      // have happened after the state dump.
      double eventPeriod = pTabDumpEvent->getPeriod();
      double now = rMolzer.eventQ.getSimTime();

      double nextTime = std::ceil(now / eventPeriod) * eventPeriod;
      
      rMolzer.eventQ.scheduleEvent(pTabDumpEvent,
				   nextTime);
    }
  };

  void
  mzrUnit::prepareToContinue(xmlpp::Element* pRootElt,
			     xmlpp::Element* pModelElt,
			     xmlpp::Element* pStreamsElt,
			     xmlpp::Element* pEventsElt,
			     std::map<std::string, std::string>& rTagToName,
			     xmlpp::Element* pTaggedSpeciesElement)
    throw(std::exception)
  {
    // Hve each tabDumpEvent schedule itself at the current simulation
    // time, which is now the time at which the state dump happened.
    std::for_each(tabDumpEvents.begin(),
		  tabDumpEvents.end(),
		  continueTabDumpEvent(rMolzer));
  }
}
