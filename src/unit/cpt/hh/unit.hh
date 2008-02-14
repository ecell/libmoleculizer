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

#ifndef CPT_UNIT_H
#define CPT_UNIT_H


#include <set>
#include <string>
#include "utl/dom.hh"
#include "cpt/inputCapabilities.hh"

namespace cpt
{
  class cptApp;

  class unit
  {
  public:

    // For debugging only, at present.  Possibly useful for other things?
    std::string name;

    cptApp& rCptApp;

    unit(const std::string& rName,
	 cptApp& refCptApp) :
      name(rName),
      rCptApp(refCptApp)
    {}

    virtual ~unit(void)
    {}

    // Keeps track of the elements that are handled by this unit,
    // so that elements don't "fall through the cracks."
    inputCapabilities inputCap;

    // The actual input parsing done by this unit.
    virtual void
    parseDomInput(xmlpp::Element* pRootElt,
		  xmlpp::Element* pModelElt,
		  xmlpp::Element* pStreamsElt,
		  xmlpp::Element* pEventsElt) throw(std::exception) = 0;

    // This allows finalization steps after parsing but before running.
    // In this phase, plexUnit runs notifications on parsed (but never created)
    // plexSpecies
    // after all parsing is complete, so that these already-defined
    // plexSpecies can be safely created just by using update.  (This is
    // still due to the clunkiness of update.)
    virtual void
    prepareToRun(xmlpp::Element* pRootElt,
		 xmlpp::Element* pModelElt,
		 xmlpp::Element* pStreamsElt,
		 xmlpp::Element* pEventsElt) throw(std::exception)
    {
      // The plex unit does a lot; the mzrUnit does a little.  So far,
      // nobody else does anything.
    }

    // This is the after-phase of parsing, analogous in parametrizer to
    // prepareToRun's role in moleculizer.
    virtual void
    prepareToDump(xmlpp::Element* pRootElt,
		  xmlpp::Element* pModelElt,
		  xmlpp::Element* pStreamsElt,
		  xmlpp::Element* pEventsElt,
		  xmlpp::Element* pTaggedSpeciesElement) throw(std::exception)
    {
      // Default implementation is to ignore the tagged-species (which most
      // units don't have) and to go on as if preparing to do a moleculizer
      // simulation.
      prepareToRun(pRootElt,
		   pModelElt,
		   pStreamsElt,
		   pEventsElt);
    }

    // This is the after-phase of parsing, analogous in continuator to
    // prepareToRun's role in moleculizer and that of prepareToDump in
    // parametrizer.
    virtual void
    prepareToContinue(xmlpp::Element* pRootElt,
		      xmlpp::Element* pModelElt,
		      xmlpp::Element* pStreamsElt,
		      xmlpp::Element* pEventsElt,
		      std::map<std::string, std::string>& rTagToName,
		      xmlpp::Element* pTaggedSpeciesElement)
      throw(std::exception)
    {
      // Default implementation is to ignore the tagged-species (which most
      // units don't have) and to go on as if preparing to do a moleculizer
      // simulation.
      prepareToRun(pRootElt,
		   pModelElt,
		   pStreamsElt,
		   pEventsElt);
    }

    // How a unit contributes its part of a state dump.
    virtual void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception) = 0;
  };
}

#endif // CPT_UNIT_H
