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

#include <set>
#include <algorithm>
#include <fstream>
#include <ctime>
#include "platform.hh"
#include "sampleDist/sampleDist.hh"
#include "domUtils/domUtils.hh"
#include "mzr/unit.hh"
#include "mzr/species.hh"
#include "mzr/reaction.hh"
#include "mzr/reactionFamily.hh"
#include "mzr/event.hh"
#include "mzr/moleculizer.hh"
#include "mzr/unitsMgr.hh"
#include "mzr/mzrUnitParse.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/linearHash.hh"

namespace mzr
{
  // Some predicates for testing that every node is within the capabilities
  // of this version of moleculizer.
  namespace
  {
    class modelNodeNotInCap :
      public std::unary_function<xmlpp::Node*, bool>
    {
      inputCapabilities& rInputCap;
    public:
      modelNodeNotInCap(inputCapabilities& rInputCapabilities) :
	rInputCap(rInputCapabilities)
      {}
    
      bool
      operator()(xmlpp::Node* pNode) const throw(mzrXcpt)
      {
	xmlpp::Element* pElt = domUtils::mustBeElementPtr(pNode);
	return (! rInputCap.handlesModelContentElt(pElt));
      }
    };

    class explicitSpeciesNodeNotInCap :
      public std::unary_function<xmlpp::Node*, bool>
    {
      inputCapabilities& rInputCap;
    public:
      explicitSpeciesNodeNotInCap(inputCapabilities& rInputCapabilities) :
	rInputCap(rInputCapabilities)
      {}
    
      bool
      operator()(xmlpp::Node* pNode) const throw(mzrXcpt)
      {
	xmlpp::Element* pElt = domUtils::mustBeElementPtr(pNode);
	return (! rInputCap.handlesExplictSpeciesContent(pElt));
      }
    };

    class speciesStreamNodeNotInCap :
      public std::unary_function<xmlpp::Node*, bool>
    {
      inputCapabilities& rInputCap;
    public:
      speciesStreamNodeNotInCap(inputCapabilities& rInputCapabilities) :
	rInputCap(rInputCapabilities)
      {}
    
      bool
      operator()(xmlpp::Node* pNode) const throw(mzrXcpt)
      {
	xmlpp::Element* pElt = domUtils::mustBeElementPtr(pNode);
	return (! rInputCap.handlesSpeciesStreamsContent(pElt));
      }
    };

    class eventNodeNotInCap :
      public std::unary_function<xmlpp::Node*, bool>
    {
      inputCapabilities& rInputCap;
    public:
      eventNodeNotInCap(inputCapabilities& rInputCapabilities) :
	rInputCap(rInputCapabilities)
      {}
    
      bool
      operator()(xmlpp::Node* pNode) const throw(mzrXcpt)
      {
	xmlpp::Element* pElt = domUtils::mustBeElementPtr(pNode);
	return (! rInputCap.handlesEventsContent(pElt));
      }
    };

    class reactionGenNotInCap :
      public std::unary_function<xmlpp::Node*, bool>
    {
      inputCapabilities& rInputCap;
    public:
      reactionGenNotInCap(inputCapabilities& rInputCapabilities) :
	rInputCap(rInputCapabilities)
      {}
    
      bool
      operator()(xmlpp::Node* pNode) const throw(mzrXcpt)
      {
	xmlpp::Element* pElt = domUtils::mustBeElementPtr(pNode);
	return (! rInputCap.handlesReactionGensContent(pElt));
      }
    };


    // For getting each unit to do its part of parsing input.
    class unitParseDomInput :
      public std::unary_function<unit*, void>
    {
      xmlpp::Element* pRootElt;
      xmlpp::Element* pModelElt;
      xmlpp::Element* pStreamsElt;
      xmlpp::Element* pEventsElt;
    public:
      unitParseDomInput(xmlpp::Element* pRootElement,
			xmlpp::Element* pModelElement,
			xmlpp::Element* pStreamsElement,
			xmlpp::Element* pEventsElement) :
	pRootElt(pRootElement),
	pModelElt(pModelElement),
	pStreamsElt(pStreamsElement),
	pEventsElt(pEventsElement)
      {}
      void
      operator()(unit* pUnit) const throw(std::exception)
      {
	pUnit->parseDomInput(pRootElt,
			     pModelElt,
			     pStreamsElt,
			     pEventsElt);
      }
    };
  }

  class moleculizer::unhandledModelContentXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingModelContentNode)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingModelContentNode)
		<< "No unit claims to handle this "
		<< pOffendingModelContentNode->get_name()
		<< " node in the model section.";
      return msgStream.str();
    }
  public:
    unhandledModelContentXcpt(xmlpp::Node* pOffendingModelContentNode) :
      mzrXcpt(mkMsg(pOffendingModelContentNode))
    {}
  };

  class moleculizer::unhandledExplicitSpeciesContentXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pBadExplicitSpeciesContentNode)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pBadExplicitSpeciesContentNode)
		<< "No unit claims to handle this "
		<< pBadExplicitSpeciesContentNode->get_name()
		<< " node in the explicit species section.";
      return msgStream.str();
    }
  public:
    unhandledExplicitSpeciesContentXcpt(xmlpp::Node* pBadExplicitSpeciesContentNode) :
      mzrXcpt(mkMsg(pBadExplicitSpeciesContentNode))
    {}
  };
    
  class moleculizer::unhandledSpeciesStreamsContentXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pBadSpeciesStreamsContentNode)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pBadSpeciesStreamsContentNode)
		<< "No unit claims to handle this "
		<< pBadSpeciesStreamsContentNode->get_name()
		<< " node in the species streams section.";
      return msgStream.str();
    }
  public:
    unhandledSpeciesStreamsContentXcpt(xmlpp::Node* pBadSpeciesStreamsContentNode) :
      mzrXcpt(mkMsg(pBadSpeciesStreamsContentNode))
    {}
  };

  class moleculizer::unhandledEventsContentXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pOffendingEventsContentNode)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingEventsContentNode)
		<< "No unit claims to handle this "
		<< pOffendingEventsContentNode->get_name()
		<< " node in the events section.";
      return msgStream.str();
    }

  public:
    unhandledEventsContentXcpt(xmlpp::Node* pOffendingEventsContentNode) :
      mzrXcpt(mkMsg(pOffendingEventsContentNode))
    {}
  };

  class moleculizer::unhandledReactionGenXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(xmlpp::Node* pBadReactionGenNode)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pBadReactionGenNode)
		<< "No unit claims to handle this "
		<< pBadReactionGenNode->get_name()
		<< " node in the reaction generators section.";
      return msgStream.str();
    }
  public:
    unhandledReactionGenXcpt(xmlpp::Node* pBadReactionGenNode) :
      mzrXcpt(mkMsg(pBadReactionGenNode))
    {
    }
  };

  void
  moleculizer::constructorPrelude(void)
  {
    // Add up the input capabilities of the units.
    pUserUnits->unionInputCaps(inputCap);

    // Formerly initialized random seed here.
  }

  void
  moleculizer::constructorCore(xmlpp::Element* pRootElement,
			       xmlpp::Element* pModelElement,
			       xmlpp::Element* pStreamsElement,
			       xmlpp::Element* pEventsElement)
    throw(std::exception)
  {
    // Verify that every element in the model section is handled by some
    // unit or another.
    xmlpp::Node::NodeList modelContentNodes
      = pModelElement->get_children();
    xmlpp::Node::NodeList::iterator iUnhandledModelContent
      = std::find_if(modelContentNodes.begin(),
		     modelContentNodes.end(),
		     modelNodeNotInCap(inputCap));
    if(modelContentNodes.end() != iUnhandledModelContent)
      throw unhandledModelContentXcpt(*iUnhandledModelContent);

    // Get the reaction-gens node, which contains kinds of reaction generators
    // introduced by units.
    xmlpp::Element* pReactionGensElt
      = domUtils::mustGetUniqueChild(pModelElement,
				     eltName::reactionGens);

    // Verify that every kind of reaction generator is handled by some
    // unit or another.
    xmlpp::Node::NodeList reactionGensContentNodes
      = pReactionGensElt->get_children();
    xmlpp::Node::NodeList::iterator iUnhandledReactionGenContent
      = std::find_if(reactionGensContentNodes.begin(),
		     reactionGensContentNodes.end(),
		     reactionGenNotInCap(inputCap));
    if(reactionGensContentNodes.end() !=
       iUnhandledReactionGenContent)
      throw unhandledReactionGenXcpt(*iUnhandledReactionGenContent);

    // Get the explicit species model node, which can contain kinds of explicit
    // species introduced by units.
    xmlpp::Element* pExplicitSpeciesElt
      = domUtils::mustGetUniqueChild(pModelElement,
				     eltName::explicitSpecies);

    // Verify that every kind of explicit species is handled by some
    // unit or another.
    xmlpp::Node::NodeList explicitSpeciesContentNodes
      = pExplicitSpeciesElt->get_children();
    xmlpp::Node::NodeList::iterator iUnhandledExplicitSpeciesContent
      = std::find_if(explicitSpeciesContentNodes.begin(),
		     explicitSpeciesContentNodes.end(),
		     explicitSpeciesNodeNotInCap(inputCap));
    if(explicitSpeciesContentNodes.end() != iUnhandledExplicitSpeciesContent)
      throw
	unhandledExplicitSpeciesContentXcpt(*iUnhandledExplicitSpeciesContent);

    // Get the speciesStreams node, which can contain kinds of species streams
    // introduced by units.
    xmlpp::Element* pSpeciesStreamsElt
      = domUtils::mustGetUniqueChild(pStreamsElement,
				     eltName::speciesStreams);
  
    // Verify that every kind of species stream is handled by some
    // unit or another.
    xmlpp::Node::NodeList speciesStreamsContentNodes
      = pSpeciesStreamsElt->get_children();
    xmlpp::Node::NodeList::iterator iUnhandledSpeciesStreamsContent
      = std::find_if(speciesStreamsContentNodes.begin(),
		     speciesStreamsContentNodes.end(),
		     speciesStreamNodeNotInCap(inputCap));
    if(speciesStreamsContentNodes.end() != iUnhandledSpeciesStreamsContent)
      throw unhandledSpeciesStreamsContentXcpt(*iUnhandledSpeciesStreamsContent);
  
    // Verify that every element in the events section is handled by some
    // unit or another.
    xmlpp::Node::NodeList eventsContentNodes
      = pEventsElement->get_children();
    xmlpp::Node::NodeList::iterator iUnhandledEventsContent
      = std::find_if(eventsContentNodes.begin(),
		     eventsContentNodes.end(),
		     eventNodeNotInCap(inputCap));
    if(eventsContentNodes.end() != iUnhandledEventsContent)
      throw unhandledEventsContentXcpt(*iUnhandledEventsContent);

    // Have each unit do its parsing thing.
    std::for_each(pUserUnits->begin(),
		  pUserUnits->end(),
		  unitParseDomInput(pRootElement,
				    pModelElement,
				    pStreamsElement,
				    pEventsElement));
  }

  moleculizer::moleculizer(void)
    throw(std::exception) :
    pUserUnits(new unitsMgr(*this))
  {}

  moleculizer::moleculizer(int argc,
			   char* argv[],
			   xmlpp::Document* pDoc)
    throw(std::exception) :
    pUserUnits(new unitsMgr(*this))
  {
    // This also initializes the random seed.
    processCommandLineArgs(argc, argv);

    // Now just does the "input capabilities" thing.
    constructorPrelude();

    // Get the basic framework of the document.
    xmlpp::Element* pRootElement
      = pDoc->get_root_node();

    // Get the unique model element.
    xmlpp::Element* pModelElement
      = domUtils::mustGetUniqueChild(pRootElement,
				     eltName::model);

    // Get the unique streams element.
    xmlpp::Element* pStreamsElement
      = domUtils::mustGetUniqueChild(pRootElement,
				     eltName::streams);

    // Get the unique events element.
    xmlpp::Element* pEventsElement
      = domUtils::mustGetUniqueChild(pRootElement,
				     eltName::events);

    // Extract model info.
    constructorCore(pRootElement,
		    pModelElement,
		    pStreamsElement,
		    pEventsElement);
  
    // Have each unit do its prepareToRun thing.
    std::for_each(pUserUnits->begin(),
		  pUserUnits->end(),
		  prepareUnitToRun(pRootElement,
				   pModelElement,
				   pStreamsElement,
				   pEventsElement));
  }

  moleculizer::~moleculizer(void)
  {
    delete pUserUnits;
  }

  void
  moleculizer::processCommandLineArgs(int argc,
				      char* argv[])
  {
    // Since this routine CAN initialize the random number generator,
    // I'm going to make it THE routine that initializes the rng.
    // It does it twice when a random seed is given on the command line.
    sampleDist::seedUniformSampler(42);

    // Skip the command name.
    argc--;
    argv++;

    // Peel off arguments one by one.
    while(0 < argc)
      {
	std::string arg(*argv);
	argv++;
	argc--;

	// The timeout time.
	if(arg == "-t")
	  {
	    if(argc <= 0) throw insufficientArgsXcpt();
	    
	    std::string timeoutString(*argv);
	    argv++;
	    argc--;

	    int timeoutSeconds = 0;
	    if(! (domUtils::stringIsInt(timeoutString,
					timeoutSeconds)
		  && (0 < timeoutSeconds)))
	      throw badPosIntArgXcpt(timeoutString);

	    pUserUnits->theMzrUnit.setTimeout((time_t) timeoutSeconds);
	  }

	// A random seed string.
	else if(arg == "-s")
	  {
	    if(argc <= 0) throw insufficientArgsXcpt();

	    std::string seedString(*argv);
	    argv++;
	    argc--;

	    linearHash lh;
	    sampleDist::seedUniformSampler(lh(seedString));
	  }

	// Turns off reaction generation, except for what must happen
	// before the simulation starts.
	else if (arg == "-g")
	  {
	    pUserUnits->theMzrUnit.generateOption = false;
	  }
	    
	else throw unrecognizedArgXcpt(arg);
      }
  }

  class unitInsertStateElements :
    public std::unary_function<unit*, void>
  {
    xmlpp::Element* pModelElt;
  public:
    unitInsertStateElements(xmlpp::Element* pModelElement) :
      pModelElt(pModelElement)
    {}

    void
    operator()(unit* pUnit) const throw(std::exception)
    {
      pUnit->insertStateElts(pModelElt);
    }
  };

  int
  moleculizer::run(void) throw(std::exception)
  {
    eventQ.run(*this,
	       pUserUnits->theMzrUnit);
    return 0;
  }

  // State dump, invoked in a dump-state event.
  xmlpp::Document*
  moleculizer::makeDomOutput(void) throw(std::exception)
  {
    xmlpp::Document* pDoc = new xmlpp::Document();

    // Create the moleculizer-state node.
    xmlpp::Element* pRootElt
      = pDoc->create_root_node(eltName::moleculizerState);

    // Add some mandatory elements.
    xmlpp::Element* pModelElt
      = pRootElt->add_child(eltName::model);
    pModelElt->add_child(eltName::unitsStates);
    pModelElt->add_child(eltName::explicitSpeciesTags);
    pModelElt->add_child(eltName::taggedSpecies);
    pModelElt->add_child(eltName::tagReactions);

    // The only reason for inserting these elements here seems
    // to have been (thought to be?) that various modules would insert
    // elements under them.
    pModelElt->add_child(eltName::time);
    pRootElt->add_child(eltName::streams);

    // Run through the units, letting each make its complete state
    // contribution, for now.
    std::for_each(pUserUnits->begin(),
		  pUserUnits->end(),
		  unitInsertStateElements(pRootElt));

    return pDoc;
  }

  // There is no .cc file corresponding to species.hh
  int species::speciesCount = 0;
}
