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

#include "utl/dom.hh"
#include "mzr/unit.hh"
#include "mzr/unitsMgr.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/parametrizer.hh"

namespace mzr
{
  class prepareUnitToDump :
    public std::unary_function<unit*, void>
  {
    xmlpp::Element* pRootElt;
    xmlpp::Element* pModelElt;
    xmlpp::Element* pStreamsElt;
    xmlpp::Element* pEventsElt;
    xmlpp::Element* pTaggedSpeciesElt;
    
  public:
    prepareUnitToDump(xmlpp::Element* pRootElement,
		      xmlpp::Element* pModelElement,
		      xmlpp::Element* pStreamsElement,
		      xmlpp::Element* pEventsElement,
		      xmlpp::Element* pTaggedSpeciesElement) :
      pRootElt(pRootElement),
      pModelElt(pModelElement),
      pStreamsElt(pStreamsElement),
      pEventsElt(pEventsElement),
      pTaggedSpeciesElt(pTaggedSpeciesElement)
    {}

    void
    operator()(unit* pUnit) const
      throw(std::exception)
    {
      pUnit->prepareToDump(pRootElt,
			   pModelElt,
			   pStreamsElt,
			   pEventsElt,
			   pTaggedSpeciesElt);
    }
  };
  
  parametrizer::parametrizer(int argc,
			     char** argv,
			     xmlpp::Document* pMoleculizerInput,
			     xmlpp::Document* pMoleculizerState)
    throw(std::exception)
  {
    // Must invoke; this initializes random seed.
    processCommandLineArgs(argc, argv);
    
    // Do the "input capabilities" thing.
    constructorPrelude();

    // Get the basic framework of moleculizer-input.
    xmlpp::Element* pInputRootElement
      = pMoleculizerInput->get_root_node();

    xmlpp::Element* pInputModelElement
      = utl::dom::mustGetUniqueChild(pInputRootElement,
				     eltName::model);
    xmlpp::Element* pInputStreamsElement
      = utl::dom::mustGetUniqueChild(pInputRootElement,
				     eltName::streams);
    xmlpp::Element* pInputEventsElement
      = utl::dom::mustGetUniqueChild(pInputRootElement,
				     eltName::events);

    // Extract model info.
    constructorCore(pInputRootElement,
		    pInputModelElement,
		    pInputStreamsElement,
		    pInputEventsElement);

    // Similar digestion of moleculizer-state.
    xmlpp::Element* pStateRootElement
      = pMoleculizerState->get_root_node();

    xmlpp::Element* pStateModelElement
      = utl::dom::mustGetUniqueChild(pStateRootElement,
				     eltName::model);
    xmlpp::Element* pStateTaggedSpeciesElement
      = utl::dom::mustGetUniqueChild(pStateModelElement,
				     eltName::taggedSpecies);

    // Have each unit do its (rather inapproriately-named) prepareToDump thing,
    // analogous to prepare-to-run in moleculizer.  The default implementation
    // is prepare-to-run.
    std::for_each(pUserUnits->begin(),
		  pUserUnits->end(),
		  prepareUnitToDump(pInputRootElement,
				    pInputModelElement,
				    pInputStreamsElement,
				    pInputEventsElement,
				    pStateTaggedSpeciesElement));
  }

  parametrizer::~parametrizer(void)
  {}

  int
  parametrizer::run(void) throw(std::exception)
  {
    xmlpp::Document* pOutputDoc = makeDomOutput();

    // This is hideous, but write_to_stream seems not to work; it produces a
    // huge, bad output file with many copies of the data in it.
    std::string output = pOutputDoc->write_to_string();

    delete pOutputDoc;

    std::cout << output;

    return 0;
  }
}
