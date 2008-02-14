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

#ifndef EVENTQUEUE_H
#define EVENTQUEUE_H

#include <map>
#include "utl/xcpt.hh"

namespace mzr
{
  class mzrEvent;
  class moleculizer;

  class eventQueue :
    public std::multimap<double, mzrEvent*>
  {
    double now;

  public:
    //! Time for an event that never happens.
    static const double never;

    eventQueue(void) :
      now(0.0)
    {}

    // To set the initial simulation time.
    void
    setSimTime(double simTime)
    {
      now = simTime;
    }

    //! Gets the current simulation time.
    double getSimTime(void)
    {
      return now;
    }
    
    //! Deschedules an event, if it is scheduled.
    void
    descheduleEvent(mzrEvent* pEvent);

    //! Schedule or reschedule an event.
    void
    scheduleEvent(mzrEvent* pEvent,
		  double time);

    /*! \brief The main simulation loop.

    The next event is picked off the end of the queue, and the current
    simulation time is set to the event's time.  Then the event's
    doEvent method is called.

    If doEvent returns event::eventResult::go, then eventQueue::run
    continues processing the next event in the queue.

    If doEvent returns event::eventResult::stop, then the simulation
    is stopped and the simulation moves on to next command (following
    the runCmd) if any.  If there are no more commands, moleculizer is
    done.  */
    void run(moleculizer& rMolzer)
      throw(utl::xcpt);
  };
}

#endif
