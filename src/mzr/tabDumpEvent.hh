/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008 The Molecular Sciences Institute.
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Contact information:
//   Larry Lok, Research Fellow          Voice: 510-981-8740
//   The Molecular Sciences Institute      Fax: 510-647-0699
//   2168 Shattuck Ave.                  Email: lok@molsci.org
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#ifndef TABDUMPEVENT_H
#define TABDUMPEVENT_H

#include "mzr/mzrEvent.hh"
#include "fnd/dumpStream.hh"
#include "mzr/mzrStream.hh"

namespace mzr
{
  // This event class is separated out because of its special
  // role in state dump and in mzrUnit.
  class tabDumpEvent :
    public mzrEvent,
    public fnd::dumpStream,
    public mzrDumpStream
  {
    // This is a self-rescheduling, periodic event.  Maybe periodic
    // should be a sub-class of event; it seems to be a recurrent enough
    // theme.
    double period;
  
  public:
    tabDumpEvent(double dumpPeriod,
		 const std::string& fileName);

    // Now use push_back to add an ordinary dumpable, and use
    // push_back in combination with addSpeciesStream to add a species
    // dumpable.

    double
    getPeriod(void)
    {
      return period;
    }

    fnd::eventResult
    happen(moleculizer& rMolzer)
      throw(std::exception);

    void
    insertTaggedDumpStreamElts(xmlpp::Element* pParent) const
      throw(std::exception);
  };
}


#endif // TABDUMPEVENT_H
