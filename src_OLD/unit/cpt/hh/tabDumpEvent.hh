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

#ifndef CPT_TABDUMPEVENT_H
#define CPT_TABDUMPEVENT_H

#include "cpt/cptEvent.hh"

#include "fnd/dumpStream.hh"
#include "cpt/stream.hh"

namespace cpt
{
  // This event class is separated out because of its special
  // role in state dump and in mzrUnit.
  class tabDumpEvent :
    public cptEvent,
    public fnd::dumpStream,
    public cptDumpStream
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
    happen(cptApp& rCptApp)
      throw(std::exception);

    void
    insertTaggedDumpStreamElts(xmlpp::Element* pParent) const
      throw(std::exception);
  };
}


#endif // CPT_TABDUMPEVENT_H
