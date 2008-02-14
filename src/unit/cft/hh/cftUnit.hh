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

#ifndef CFT_CFTUNIT_H
#define CFT_CFTUNIT_H

#include "cpt/cptUnit.hh"
#include "cml/cmlUnit.hh"
#include "clx/clxUnit.hh"
#include "cft/cftEltName.hh"

namespace cft
{
  class cftUnit : 
    public cpt::unit
  {
  public:
    cpt::cptUnit& rCptUnit;
    cml::cmlUnit& rCmlUnit;
    clx::clxUnit& rClxUnit;

    cftUnit(cpt::cptApp& rCptApp,
	    cpt::cptUnit& refCptUnit,
	    cml::cmlUnit& refCmlUnit,
	    clx::clxUnit& refClxUnit) :
      cpt::unit("cft",
		rCptApp),
      rCptUnit(refCptUnit),
      rCmlUnit(refCmlUnit),
      rClxUnit(refClxUnit)
    {
      // Register reaction generator names.  This is used by the parser
      // to verify that all elements appearing in the input file are "claimed"
      // by some module.  Why did I think that was a good idea?  Now, "open
      // schema" i.e. pay no attention to extraneous matter while parsing.
      // seems more rational to me.
      inputCap.addReactionGenName(eltName::omniGen);
      inputCap.addReactionGenName(eltName::uniMolGen);

      // Register the enabling complexes for omniRxnGen generators
      // as omniplexes, for processing by the clxUnit.
      const std::string slash("/");
      std::ostringstream omniGensXpath;
      omniGensXpath << cpt::eltName::model
		    << slash
		    << cpt::eltName::reactionGens
		    << slash
		    << eltName::omniGen
		    << slash
		    << eltName::enablingOmniplex;
      rClxUnit.addOmniXpath(omniGensXpath.str());
    }

    // The input parsing routine for this unit.
    void
    parseDomInput(xmlpp::Element* pRootElement,
		  xmlpp::Element* pModelElement,
		  xmlpp::Element* pStreamsElement,
		  xmlpp::Element* pEventsElement)
      throw(std::exception);

    // The state output routine for this unit.
    void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
    {
      // For now, this unit just provides reaction generators, which
      // for now I'm not dumping in "state dump."
    }
  };
}

#endif // CFT_CFTUNIT_H
