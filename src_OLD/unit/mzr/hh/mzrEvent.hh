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

#ifndef MZR_MZREVENT_H
#define MZR_MZREVENT_H

#include "fnd/event.hh"
#include "mzr/eventQueue.hh"

namespace mzr
{
  class moleculizer;
  
  // Note that continuator and parametrizer inherit from moleculizer.
  class mzrEvent :
    public fnd::event<moleculizer>
  {
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

  public:
    mzrEvent(void) :
      scheduled(false)
    {}

    virtual
    ~mzrEvent(void)
    {}

    void
    schedule(eventQueue::iterator iEntry)
    {
      iQueueEntry = iEntry;
      scheduled = true;
    }

    // Returns true if the event was scheduled, and returns the
    // event queue entry at which it was scheduled in that case.
    bool
    deschedule(eventQueue::iterator& riEntry)
    {
      if(scheduled)
	{
	  scheduled = false;
	  riEntry = iQueueEntry;
	  return true;
	}
      else return false;
    }

    /*! Uses accessors to eventQueue entry to get current scheduled time.
      If the event is not in the queue, then event::never is returned. */
    double
    getTime(void)
    {
      return scheduled ? iQueueEntry->first : eventQueue::never;
    }
  };
}

#endif // MZR_MZREVENT_H
