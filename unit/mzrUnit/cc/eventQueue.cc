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
#include "mzr/eventQueue.hh"
#include "mzr/event.hh"
#include "mzr/moleculizer.hh"
#include "mzr/mzrUnit.hh"

namespace mzr
{
  class eventInPastXcpt : public mzrXcpt
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
      mzrXcpt(mkMsg(curTime,
		    badEventTime))
    {}
  };
  
  
  // Deschedule an event.
  void
  eventQueue::descheduleEvent(event* pEvent)
  {
    if(pEvent->scheduled)
      {
	erase(pEvent->iQueueEntry);
	pEvent->scheduled = false;
      }
  }

  // Schedule or reschedule an event.
  void
  eventQueue::scheduleEvent(event* pEvent,
			    double time)
  {
    // Maybe I should test that time >= now?

    // Deschedule the event if it is scheduled.
    descheduleEvent(pEvent);

    // Insert the event into the schedule.
    pEvent->iQueueEntry  = insert(std::make_pair(time, pEvent));
    pEvent->scheduled = true;
  }

  void
  eventQueue::run(moleculizer& rMolzer,
		  mzrUnit& rMzrUnit)
  {
    event::eventResult result = event::go;
    while(result == event::go)
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
	    if(event::never == now)
	      throw nowIsNeverXcpt();

	    // Remove the event from the event schedule.
	    // This has to be done before executing the event,
	    // because some events (dump) reschedule themselves.
	    event* pNextEvent = iTimeEventPair->second;
	    descheduleEvent(pNextEvent);

	    // Do the next event.
	    result = pNextEvent->doEvent(rMolzer,
					 rMzrUnit);
	  }
      }
  }
}
