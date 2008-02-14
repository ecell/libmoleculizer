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

#include <time.h>
#include <limits>
#include "mzr/eventQueue.hh"
#include "mzr/mzrEvent.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/unitsMgr.hh"
#include "mzr/moleculizer.hh"

namespace mzr
{
  class timeLimitExpiredXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(double clockTimeLimit,
	  double curSimTime)
    {
      std::ostringstream msgStream;
      msgStream << "Clock time limit of "
		<< clockTimeLimit
		<< " seconds exceeded at simulation time "
		<< curSimTime
		<< ".";
      return msgStream.str();
    }
  public:
    timeLimitExpiredXcpt(double clockTimeLimit,
			 double curSimTime) :
      utl::xcpt(mkMsg(clockTimeLimit,
		      curSimTime))
    {}
  };

  class eventQueueEmptyXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(double curSimTime)
    {
      std::ostringstream msgStream;
      msgStream << "Event queue empty at simulation time "
		<< curSimTime
		<< ".";
      return msgStream.str();
    }
  public:
    eventQueueEmptyXcpt(double curSimTime) :
      utl::xcpt(mkMsg(curSimTime))
    {}
  };

  class nowIsNeverXcpt :
    public utl::xcpt
  {
  public:
    nowIsNeverXcpt(void) :
      utl::xcpt("Earliest event is at 'infinite' simulation time.")
    {}
  };

  class eventInPastXcpt :
    public utl::xcpt
  {
    static std::string
    mkMsg(double curTime,
	  double badEventTime)
    {
      std::ostringstream msgStream;
      msgStream << "At simulation time "
		<< curTime
		<< ", next event is at past time "
		<< badEventTime
		<< ".";
      return msgStream.str();
    }
  public:
    eventInPastXcpt(double curTime,
		    double badEventTime) :
      utl::xcpt(mkMsg(curTime,
		      badEventTime))
    {}
  };
  
  // Deschedule an event.
  void
  eventQueue::descheduleEvent(mzrEvent* pEvent)
  {
    iterator iEvent;
    if(pEvent->deschedule(iEvent)) erase(iEvent);
  }

  // Schedule or reschedule an event.
  void
  eventQueue::scheduleEvent(mzrEvent* pEvent,
			    double time)
  {
    // Deschedule the event if it is scheduled.
    descheduleEvent(pEvent);

    // Insert the event into the schedule, and update the
    // event's internal queue iterator.
    pEvent->schedule(insert(std::make_pair(time, pEvent)));
  }

  void
  eventQueue::run(moleculizer& rMolzer)
    throw(utl::xcpt)
  {
    mzrUnit& rMzrUnit = *(rMolzer.pUserUnits->pMzrUnit);
    
    fnd::eventResult result = fnd::go;
    while(result == fnd::go)
      {
	// Check the elapsed time, if called for, to see if the clock
	// time limit has expired.
	if(0 < rMzrUnit.secondsLimit)
	  {
	    time_t secondsNow = time(0);
	    if(rMzrUnit.secondsLimit + rMzrUnit.startSeconds < secondsNow)
	      {
		// Write the state (to a test-file for now)
		xmlpp::Document* pStateDoc = rMolzer.makeDomOutput();
		pStateDoc->write_to_file_formatted("timeout-state.xml");
		delete pStateDoc;

		// End the simulation.
		throw timeLimitExpiredXcpt(rMzrUnit.secondsLimit,
					   now);
	      }
	  }

	// Get iterator to first time/event pair in the event
	// schedule.
	iterator iTimeEventPair = begin();

	// Is there really an event to execute?
	if(end() == iTimeEventPair)
	  throw eventQueueEmptyXcpt(now);

	double nextNow = iTimeEventPair->first;

	if(nextNow < now)
	  {
	    throw eventInPastXcpt(now,
				  nextNow);
	  }
	else
	  {
	    // Reset "now" to the time of the next event.
	    now = iTimeEventPair->first;

	    // Check that "now" is not never.
	    if(never == now)
	      throw nowIsNeverXcpt();

	    // Remove the event from the event schedule.
	    // This has to be done before executing the event,
	    // because some events (dump) reschedule themselves.
	    mzrEvent* pNextEvent = iTimeEventPair->second;
	    descheduleEvent(pNextEvent);

	    // Do the next event.
	    result = pNextEvent->happen(rMolzer);
	  }
      }
  }

  const double eventQueue::never = std::numeric_limits<double>::max();
}
