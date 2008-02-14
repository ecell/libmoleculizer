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

#ifndef MZR_CREATEEVENT_H
#define MZR_CREATEEVENT_H

#include "mzr/mzrEvent.hh"

namespace mzr
{
  class mzrSpecies;
  class mzrUnit;
  
  /*! \ingroup eventGroup
    \brief Abstract base class for all events.

    Basically, everything that happens in the simulation is an event.

    Aside from the pure virtual doEvent method, which can do just about
    anything, the event is enabled to keep track of its entry in the
    global queue of events.  This allows one to determine if the event
    is scheduled to happen and at what simulation time. */

  /*! \ingroup eventGroup
    \brief Event to create species members. */
  class createEvent :
    public mzrEvent
  {
    mzrSpecies* pSpeciesToCreate;
    int howMany;

    mzrUnit& rMzrUnit;

  public:
    createEvent(mzrSpecies* pSpecies,
		int count,
		mzrUnit& refMzrUnit);

    fnd::eventResult
    happen(moleculizer& rMolzer)
      throw(std::exception);
  };
}

#endif // MZR_CREATEEVENT_H
