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
#include "mzr/continuator.hh"

namespace mzr
{
  class prepareUnitToContinue :
    public std::unary_function<unit*, void>
  {
    xmlpp::Element* pRootElt;
    xmlpp::Element* pModelElt;
    xmlpp::Element* pStreamsElt;
    xmlpp::Element* pEventsElt;
    std::map<std::string, std::string>& rTagToName;
    xmlpp::Element* pTaggedSpeciesElt;
    
  public:
    prepareUnitToContinue(xmlpp::Element* pRootElement,
			  xmlpp::Element* pModelElement,
			  xmlpp::Element* pStreamsElement,
			  xmlpp::Element* pEventsElement,
			  std::map<std::string, std::string>& refTagToName,
			  xmlpp::Element* pTaggedSpeciesElement) :
      pRootElt(pRootElement),
      pModelElt(pModelElement),
      pStreamsElt(pStreamsElement),
      pEventsElt(pEventsElement),
      rTagToName(refTagToName),
      pTaggedSpeciesElt(pTaggedSpeciesElement)
    {}

    void
    operator()(unit* pUnit) const
      throw(std::exception)
    {
      pUnit->prepareToContinue(pRootElt,
			       pModelElt,
			       pStreamsElt,
			       pEventsElt,
			       rTagToName,
			       pTaggedSpeciesElt);
    }
  };
  
  class getTagName : public
  std::unary_function<xmlpp::Node*, std::pair<std::string, std::string> >
  {
  public:
    std::pair<std::string, std::string>
    operator()(xmlpp::Node* pExplicitSpeciesTagNode) const
      throw(std::exception)
    {
      xmlpp::Element* pExplicitSpeciesTagElt
	= utl::dom::mustBeElementPtr(pExplicitSpeciesTagNode);

      std::string speciesName
	= utl::dom::mustGetAttrString(pExplicitSpeciesTagElt,
				      eltName::explicitSpeciesTag_nameAttr);

      std::string speciesTag
	= utl::dom::mustGetAttrString(pExplicitSpeciesTagElt,
				      eltName::explicitSpeciesTag_tagAttr);

      return std::pair<std::string, std::string>(speciesTag, speciesName);
    }
  };
  
  continuator::continuator(int argc,
			   char** argv,
			   xmlpp::Document* pMoleculizerInput,
			   xmlpp::Document* pMoleculizerState)
    throw(std::exception)
  {
    // This also initializes random seed.
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

    // Similar digestion of moleculizer-state.
    xmlpp::Element* pStateRootElement
      = pMoleculizerState->get_root_node();

    xmlpp::Element* pStateModelElement
      = utl::dom::mustGetUniqueChild(pStateRootElement,
				     eltName::model);
    xmlpp::Element* pStateTaggedSpeciesElement
      = utl::dom::mustGetUniqueChild(pStateModelElement,
				     eltName::taggedSpecies);
    xmlpp::Element* pTimeElement
      = utl::dom::mustGetUniqueChild(pStateModelElement,
				     eltName::time);

    // Get the time at which the state was dumped.
    double dumpTime
      = utl::dom::mustGetAttrDouble(pTimeElement,
				    eltName::time_secondsAttr);

    // Set the current simulation time to the dump time.
    eventQ.setSimTime(dumpTime);

    // Extract model info.  This has to be done after the simulation time is
    // restored, so that explicit events in the moleculizer-input that have
    // already happened will not be scheduled.
    constructorCore(pInputRootElement,
		    pInputModelElement,
		    pInputStreamsElement,
		    pInputEventsElement);

    // Get the map from dump tag to explicit species names.
    // The stochastirator species are known only by name, so that
    // to match the dump tag of a stochastirator species up with
    // the species, we need this correspondence.
    xmlpp::Element* pExplicitSpeciesTagsElt
      = utl::dom::mustGetUniqueChild(pStateModelElement,
				     eltName::explicitSpeciesTags);
    xmlpp::Node::NodeList explicitSpeciesTagNodes
      = pExplicitSpeciesTagsElt->get_children(eltName::explicitSpeciesTag);
    std::map<std::string, std::string> tagToName;
    std::transform(explicitSpeciesTagNodes.begin(),
		   explicitSpeciesTagNodes.end(),
		   std::inserter(tagToName,
				 tagToName.begin()),
		   getTagName());

    // Have each unit do its "prepareToContinue" thing.
    std::for_each(pUserUnits->begin(),
		  pUserUnits->end(),
		  prepareUnitToContinue(pInputRootElement,
					pInputModelElement,
					pInputStreamsElement,
					pInputEventsElement,
					tagToName,
					pStateTaggedSpeciesElement));
  }
}
