//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                                                                          
//                                                                          
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published 
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//    
// END HEADER
// 
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//              
//

#ifndef EVENTQUEUE_H
#define EVENTQUEUE_H

#include <map>
#include "utl/xcpt.hh"
#include "mzr/mzrEvent.hh"

namespace mzr
{
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
