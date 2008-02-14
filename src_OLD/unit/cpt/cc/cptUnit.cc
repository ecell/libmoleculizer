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
#include <cmath>
#include "utl/badNNIntArgXcpt.hh"
#include "utl/unkArgXcpt.hh"
#include "utl/linearHash.hh"
#include "cpt/cptEltName.hh"
#include "cpt/cptUnit.hh"
#include "cpt/unitsMgr.hh"
#include "cpt/cptApp.hh"
#include "cpt/tabDumpEvent.hh"
#include "cpt/dupSpeciesNameXcpt.hh"
#include "cpt/unkSpeciesXcpt.hh"
#include "cpt/dupDumpableNameXcpt.hh"
#include "cpt/unkDumpableXcpt.hh"
#include "cpt/clockDumpable.hh"
#include "cpt/secondsDumpable.hh"
#include "cpt/reactEventCountDumpable.hh"
#include "cpt/reactCountDumpable.hh"
#include "cpt/speciesCountDumpable.hh"
#include "cpt/simTimeDumpable.hh"
#include "cpt/volumeDumpable.hh"

namespace cpt
{
  void
  cptUnit::
  mustAddNamedSpecies(const std::string& rSpeciesName,
		      globalSpecies* pSpecies,
		      xmlpp::Node* pRequestingNode)
    throw(utl::xcpt)
  {
    if(! addNamedSpecies(rSpeciesName,
			 pSpecies))
      throw dupSpeciesNameXcpt(rSpeciesName,
			       pRequestingNode);
  }
  
  void
  cptUnit::
  mustNameSpecies(const std::string& rSpeciesName,
		  globalSpecies* pSpecies,
		  xmlpp::Node* pRequestingNode)
    throw(utl::xcpt)
  {
    if((! speciesByName.addEntry(rSpeciesName,
				 pSpecies))
       || (! addSpeciesDumpable(new singleGlobalSpeciesDumpable(rSpeciesName,
								pSpecies))))
      {
	throw dupSpeciesNameXcpt(rSpeciesName,
				 pRequestingNode);
      }
  }

  globalSpecies*
  cptUnit::
  mustFindSpecies(const std::string& rSpeciesName,
		  xmlpp::Node* pRequestingNode) const
    throw(utl::xcpt)
  {
    globalSpecies* pSpecies = findSpecies(rSpeciesName);

    if(! pSpecies)
      throw unkSpeciesXcpt(rSpeciesName,
			   pRequestingNode);
    return pSpecies;
  }

  void
  cptUnit::
  mustAddDumpable(fnd::basicDumpable* pDumpable,
		  xmlpp::Node* pRequestingNode)
    throw(utl::xcpt)
  {
    if(! addDumpable(pDumpable))
      throw dupDumpableNameXcpt(pDumpable->getName(),
				pRequestingNode);
  }
  
  fnd::basicDumpable*
  cptUnit::
  mustFindDumpable(const std::string& rDumpableName,
		   xmlpp::Node* pRequestingNode) const
    throw(utl::xcpt)
  {
    fnd::basicDumpable* pDumpable
      = findDumpable(rDumpableName);

    if(! pDumpable)
      throw unkDumpableXcpt(rDumpableName,
			    pRequestingNode);

    return pDumpable;
  }

  cptUnit::
  cptUnit(cptApp& rCptApp) :
    unit("cpt",
	 rCptApp),
    generateOption(true),
    generateOk(true),
    generateDepth(0),
    startSeconds(time(0)),
    secondsLimit(0)
  {
    // Model elements whose contents are parsed by some unit
    // or another, as determined by moleculizer::parseDomInput.
    inputCap.addModelContentName(eltName::reactionGens);
    inputCap.addModelContentName(eltName::explicitSpecies);

    // Model elements that this unit actually processes.
    inputCap.addModelContentName(eltName::explicitReactions);
    inputCap.addModelContentName(eltName::compartmentGraph);

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
    addDumpable(new simTimeDumpable(rCptApp));
    addDumpable(new volumeDumpable(getCompartmentGraph()));
  }

  void
  cptUnit::
  parseCommandLine(int argc,
		   char* argv[])
  {
    // Skip the command name.
    argc--;
    argv++;

    // Peel off arguments one by one.
    while(0 < argc)
      {
	std::string arg(*argv);
	argc--;
	argv++;

	// Command line argument controlling the depth to probe
	// in reaction network generation.
	if(arg == "-d")
	  {
	    if(argc < 1) throw utl::insuffArgsXcpt::general();

	    arg = *argv;
	    argc--;
	    argv++;

	    if(! (utl::stringIsInt(arg,
				   generateDepth)
		  && (0 <= generateDepth)))
	      throw utl::badNNIntArgXcpt(arg);
	  }

	// Turns off reaction generation, except for what must happen
	// before the simulation starts.
	else if (arg == "-g")
	  {
	    generateOption = false;
	  }

	// A random seed string.
	else if(arg == "-s")
	  {
	    if(argc <= 0) throw utl::insuffArgsXcpt::general();

	    std::string seedString(*argv);
	    argv++;
	    argc--;

	    utl::linearHash lh;
	    getRng().seed(lh(seedString));
	  }

	else
	  {
	    throw utl::unkArgXcpt(arg);
	  }
      }
  }

  // Cause tabDumpEvent to emit its header to its file and does
  // the initial scheduling of the tabDumpEvent.
  class initScheduleTabDumpEvent :
    public std::unary_function<tabDumpEvent*, void>
  {
    cptApp& rApp;

  public:
    initScheduleTabDumpEvent(cptApp& rCptApp) :
      rApp(rCptApp)
    {}

    void
    operator()(tabDumpEvent* pTabDumpEvent) const
    {
      pTabDumpEvent->init();
      
      rApp.getQueue().scheduleEvent(pTabDumpEvent,
				    rApp.getSimTime());
    }
  };

  void
  cptUnit::prepareToRun(xmlpp::Element* pRootElt,
			xmlpp::Element* pModelElt,
			xmlpp::Element* pStreamsElt,
			xmlpp::Element* pEventsElt) throw(std::exception)
  {
    // Have each tabDumpEvent emit its header line, and schedule it for time 0,
    // after which it will reschedule itself periodically.  .
    std::for_each(tabDumpEvents.begin(),
		  tabDumpEvents.end(),
		  initScheduleTabDumpEvent(rCptApp));

    // Shut off reaction generation if called for by program option.
    generateOk = generateOption;
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
    cptApp& rApp;

  public:
    continueTabDumpEvent(cptApp& rCptApp) :
      rApp(rCptApp)
    {
    }

    void
    operator()(tabDumpEvent* pTabDumpEvent) const
    {
      pTabDumpEvent->init();
      
      // Schedule the dump event at the next time that it would
      // have happened after the state dump.
      double eventPeriod = pTabDumpEvent->getPeriod();
      double now = rApp.getSimTime();

      double nextTime = std::ceil(now / eventPeriod) * eventPeriod;
      
      rApp.getQueue().scheduleEvent(pTabDumpEvent,
				    nextTime);
    }
  };

  void
  cptUnit::prepareToContinue(xmlpp::Element* pRootElt,
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
		  continueTabDumpEvent(rCptApp));

    // Shut off reaction generation if called for by program option.
    generateOk = generateOption;
  }
}
