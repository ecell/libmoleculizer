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

#ifndef EVENT_H
#define EVENT_H

/*! \defgroup eventGroup Events and the event queue
  \ingroup mzrGroup
  \brief Events and the structure that keeps them in order. */

/*! \file event.hh
  \ingroup eventGroup
  \brief Declares event base class and several built-in descendents. */

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <float.h>
#include "mzr/eventQueue.hh"
#include "mzr/dumpable.hh"

namespace mzr
{
  class mzrUnit;
  
  /*! \ingroup eventGroup
    \brief Abstract base class for all events.

    Basically, everything that happens in the simulation is an event.

    Aside from the pure virtual doEvent method, which can do just about
    anything, the event is enabled to keep track of its entry in the
    global queue of events.  This allows one to determine if the event
    is scheduled to happen and at what simulation time. */
  class event
  {
  public:
    //! Time for an event that never happens.
    static const double never = DBL_MAX;

    /*! \brief Return value of event::doEvent.

    Tells Moleculizer whether to continue or not after the event has
    happened.  Only stopEvent tells it to stop. */
    enum eventResult
      {
	go,
	stop
      };

    /*! \name Event queue entry access.
    
    Access is through an iterator to the event's entry in the event schedule.
    The bool just tells if the iterator is valid or not.

    Using end() might have been a better way, but this way allows a
    default constructor, not needing access to the event queue to
    initialize the iterator. */
    //@{
    eventQueue::iterator iQueueEntry;
    bool scheduled;
    //@}

    event(void) :
      scheduled(false)
    {}

    virtual
    ~event(void)
    {}

    /*! This is where the event does its real work. */
    virtual eventResult
    doEvent(moleculizer& rMolzer,
	    mzrUnit& rMzrUnit) = 0;

    /*! Uses accessors to eventQueue entry to get current scheduled time.
      If the event is not in the queue, then event::never is returned. */
    double getTime(void);
  };

  /*! \ingroup eventGroup
    \brief Event to create species members. */
  class createEvent : public event
  {
    species* pSpeciesToCreate;
    int howMany;

    mzrUnit& rMzrUnit;

  public:
    createEvent(species* pSpecies,
		int count,
		mzrUnit& refMzrUnit);

    eventResult
    doEvent(moleculizer& rMolzer,
	    mzrUnit& rMzrUnit);
  };

  /*! \ingroup eventGroup
    \brief Event to change the volume once. */
  class volumeEvent : public event
  {
    mzrUnit& rMzrUnit;
    double vol;

  public:
    volumeEvent(mzrUnit& refMzrUnit,
		double targetVolume) :
      rMzrUnit(refMzrUnit),
      vol(targetVolume)
    {}
  
    virtual eventResult
    doEvent(moleculizer& rMolzer,
	    mzrUnit& rMzrUnit);
  };

  /*! \ingroup eventGroup
    \brief Periodic event for growth or shrinkage of volume. */
  class growEvent : public event
  {
    mzrUnit& rMzrUnit;
    double factor;
    double period;
  public:
    growEvent(mzrUnit& refMzrUnit,
	      double growthFactor,
	      double schedulingPeriod) :
      rMzrUnit(refMzrUnit),
      factor(growthFactor),
      period(schedulingPeriod)
    {}

    eventResult
    doEvent(moleculizer& rMolzer,
	    mzrUnit& rMzrUnit);
  };

  /*! \ingroup eventGroup
    \brief Event that stops the simulation.

    It doesn't end the program; the simulation may be restarted.
    with a runCmd.  This event is scheduled by runCmd. */
  class stopEvent : public event
  {
  public:
    eventResult
    doEvent(moleculizer& rMolzer,
	    mzrUnit& rMzrUnit);
  };

  /*! \ingroup eventGroup
    \brief Event that periodically dumps simulation data.

    This event dumps a tab-separated list of dumpables, headed by the
    simulation time and followed by a newline.  This is for dumping
    things to be plotted, put into a spreadsheet, etc. */
  class tabDumpEvent :
    public event
  {
    // The dumpables that will appear as columns in the output
    // file.
    std::vector<dumpable*> dumpables;
    // These play a special role in state dump; they're the
    // things I expect to be dumpable in other programs.
    // When a state dump happens, these get listed out.
    std::vector<speciesDumpable*> speciesDumpables;
    
    std::ostream* pOs;
    std::ofstream* pFileStream;
    // In order to make a sensible state-dump.
    std::string fileName;

    // This is a self-rescheduling, periodic event.  Maybe periodic
    // should be a sub-class of event; it seems to be a recurrent enough
    // theme.
    double period;
  
  public:
    tabDumpEvent(double dumpPeriod,
		 const std::string& fileName);

    ~tabDumpEvent(void)
    {
      delete pFileStream;
    }

    std::ostream& getOstream(void) const
    {
      return *pOs;
    }

    double
    getPeriod(void) const
    {
      return period;
    }

    void
    addDumpable(dumpable* pDumpable)
    {
      dumpables.push_back(pDumpable);
    }

    // For handling species dumpables in a way that will allow
    // them to be found when the time comes to do a state dump.
    void
    addSpeciesDumpable(speciesDumpable* pSpeciesDumpable)
    {
      addDumpable(pSpeciesDumpable);
      speciesDumpables.push_back(pSpeciesDumpable);
    }

    // Initializes output stream and writes column headers.
    void
    init(void);

    // Inserts "tagged-dump-stream" into state dump.
    void
    insertTaggedDumpStreamElts(xmlpp::Element* pParent) const
      throw(std::exception);

    eventResult
    doEvent(moleculizer& rMolzer,
	    mzrUnit& rMzrUnit);
  };

  // Dump of state for export to other applications or other processing.
  class dumpStateEvent : public event
  {
    std::string outputFileName;
  public:
    dumpStateEvent(const std::string& rOutputFileName) :
      outputFileName(rOutputFileName)
    {}
  
    eventResult
    doEvent(moleculizer& rMolzer,
	    mzrUnit& rMzrUnit) throw(std::exception);
  };

  // Event that stops reaction generation.
  class noReactEvent : public event
  {
    mzrUnit& rMzrUnit;
  public:
    noReactEvent(mzrUnit& refMzrUnit) :
      rMzrUnit(refMzrUnit)
    {}
    
    eventResult
    doEvent(moleculizer& rMolzer,
	    mzrUnit& rMzrUnit);
  };
}

#endif
