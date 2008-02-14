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

#include "mzr/mzrEltName.hh"
#include "mzr/moleculizer.hh"
#include "mzr/mzrReaction.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/stopEventInPastXcpt.hh"
#include "mzr/unkSpeciesXcpt.hh"
#include "mzr/unkStatStreamXcpt.hh"
#include "mzr/noStopEventWarning.hh"
#include "mzr/createEvent.hh"
#include "mzr/stopEvent.hh"
#include "mzr/dumpStateEvent.hh"

namespace mzr
{
  class addProductSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzrUnit& rMzrUnit;
    mzrReaction* pReaction;
  public:
    addProductSpecies(mzrUnit& refMzrUnit,
		      mzrReaction* pParsedReaction) :
      rMzrUnit(refMzrUnit),
      pReaction(pParsedReaction)
    {}

    void
    operator()(xmlpp::Node* pProductSpeciesRefNode) const throw(std::exception)
    {
      xmlpp::Element* pProductSpeciesRefElt
	= utl::dom::mustBeElementPtr(pProductSpeciesRefNode);

      // Get the name of the substrate species.
      std::string speciesName
	= utl::dom::mustGetAttrString(pProductSpeciesRefElt,
				      eltName::productSpeciesRef_nameAttr);

      // Look up the species in moleculizer's BAD BAD BAD static catalog of
      // named species.
      mzrSpecies* pSpecies
	= rMzrUnit.mustFindSpecies(speciesName,
				   pProductSpeciesRefElt);

      // Get the multiplicity of the product, which can be positive or
      // negative.
      int multiplicity
	= utl::dom::mustGetAttrInt(pProductSpeciesRefElt,
				   eltName::productSpeciesRef_multAttr);

      // Add the product species/multiplicity to the reaction.
      pReaction->addProduct(pSpecies,
			    multiplicity);
    }
  };

  class addSubstrateSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzrUnit& rMzrUnit;
    mzrReaction* pReaction;
  public:
    addSubstrateSpecies(mzrUnit& refMzrUnit,
			mzrReaction* pParsedReaction) :
      rMzrUnit(refMzrUnit),
      pReaction(pParsedReaction)
    {}

    void
    operator()(xmlpp::Node* pSubstrateSpeciesRefNode) const throw(std::exception)
    {
      // Cast the Node* to Element*, probably unnecessarily dynamically.
      xmlpp::Element* pSubstrateSpeciesRefElt
	= utl::dom::mustBeElementPtr(pSubstrateSpeciesRefNode);

      // Get the name of the substrate species.
      std::string speciesName
	= utl::dom::mustGetAttrString(pSubstrateSpeciesRefElt,
				      eltName::substrateSpeciesRef_nameAttr);

      // Look up the species in moleculizer's BAD BAD BAD static catalog of
      // named species.
      mzrSpecies* pSpecies
	= rMzrUnit.mustFindSpecies(speciesName,
				   pSubstrateSpeciesRefElt);

      // Get the multiplicity of the substrate, which must be positive.
      int multiplicity
	= utl::dom::mustGetAttrPosInt(pSubstrateSpeciesRefElt,
				      eltName::substrateSpeciesRef_multAttr);

      // Add the substrate/multiplicity to the reaction.
      pReaction->addReactant(pSpecies,
			     multiplicity);
    }
  };

  class installReaction :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzrUnit& rMzrUnit;
  public:
    installReaction(mzrUnit& refMzrUnit) :
      rMzrUnit(refMzrUnit)
    {}
    
    void
    operator()(const xmlpp::Node* pReactionNode) const throw(std::exception)
    {
      mzrReaction* pParsedReaction
	= new mzrReaction(rMzrUnit.globalVars.begin(),
			  rMzrUnit.globalVars.end());
    
      // Get the list of substrate nodes.
      xmlpp::Node::NodeList substrateSpeciesRefNodes
	= pReactionNode->get_children(eltName::substrateSpeciesRef);
      // Add the substrates to the reaction.
      std::for_each(substrateSpeciesRefNodes.begin(),
		    substrateSpeciesRefNodes.end(),
		    addSubstrateSpecies(rMzrUnit,
					pParsedReaction));

      // Get the list of product nodes.
      xmlpp::Node::NodeList productSpeciesRefNodes
	= pReactionNode->get_children(eltName::productSpeciesRef);
      // Add the products to the reaction.
      std::for_each(productSpeciesRefNodes.begin(),
		    productSpeciesRefNodes.end(),
		    addProductSpecies(rMzrUnit,
				      pParsedReaction));

      // Get the reaction rate element.
      xmlpp::Element* pRateElt
	= utl::dom::mustGetUniqueChild(pReactionNode,
				       eltName::rate);
      double rate
	= utl::dom::mustGetAttrDouble(pRateElt,
				      eltName::rate_valueAttr);
      pParsedReaction->setRate(rate);

      // Add the reaction to its destruction pit.
      rMzrUnit.addUserReaction(pParsedReaction);
    }
  };

  class resolveSpeciesStreamRef : public
  std::unary_function<xmlpp::Node*, fnd::dumpable<fnd::basicDumpable::dumpArg>*>
  {
    mzrUnit& rMzrUnit;
  public:
    resolveSpeciesStreamRef(mzrUnit& refMzrUnit) :
      rMzrUnit(refMzrUnit)
    {}

    fnd::dumpable<fnd::basicDumpable::dumpArg>*
    operator()(xmlpp::Node* pSpeciesStreamRefNode) const
      throw(std::exception)
    {
      xmlpp::Element* pSpeciesStreamRefElt
	= utl::dom::mustBeElementPtr(pSpeciesStreamRefNode);

      std::string dumpableName
	= utl::dom::mustGetAttrString(pSpeciesStreamRefElt,
				      eltName::speciesStreamRef_nameAttr);

      // Look up the speciesDumpable.
      fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable
	= rMzrUnit.mustFindDumpable(dumpableName,
				    pSpeciesStreamRefElt);

      return pDumpable;
    }
  };

  // Generates dumpables for named species.
  class resolveSpeciesRef : public
  std::unary_function<xmlpp::Node*, singleSpeciesDumpable<mzrSpecies>*>
  {
    mzrUnit& rMzrUnit;
  public:
    resolveSpeciesRef(mzrUnit& refMzrUnit) :
      rMzrUnit(refMzrUnit)
    {}

    singleSpeciesDumpable<mzrSpecies>*
    operator()(xmlpp::Node* pSpeciesRefNode) const
      throw(std::exception)
    {
      xmlpp::Element* pSpeciesRefElt
	= utl::dom::mustBeElementPtr(pSpeciesRefNode);

      std::string speciesName
	= utl::dom::mustGetAttrString(pSpeciesRefElt,
				      eltName::speciesRef_nameAttr);

      mzrSpecies* pSpecies = rMzrUnit.findSpecies(speciesName);
      if(0 == pSpecies)
	throw unkSpeciesXcpt(speciesName,
			     pSpeciesRefElt);

      // Dumpables are all memory managed by mzrUnit.
      singleSpeciesDumpable<mzrSpecies>* pDumpable 
	= new singleSpeciesDumpable<mzrSpecies>(speciesName,
						pSpecies);
      rMzrUnit.addDumpable(pDumpable);

      // Dumpables that dump populations of species or sums of populations
      // of species are special for state dump.  (To reflect this, they inherit
      // from mzrSpeciesStream.
      rMzrUnit.addSpeciesStream(pDumpable);

      return pDumpable;
    }
  };

  class resolveStatStreamRef : public
  std::unary_function<xmlpp::Node*, fnd::dumpable<fnd::basicDumpable::dumpArg>*>
  {
    mzrUnit& rMzrUnit;
  public:
    resolveStatStreamRef(mzrUnit& refMzrUnit) :
      rMzrUnit(refMzrUnit)
    {}

    fnd::dumpable<fnd::basicDumpable::dumpArg>*
    operator()(xmlpp::Node* pStatStreamRefNode) const throw(std::exception)
    {
      xmlpp::Element* pStatStreamRefElt
	= utl::dom::mustBeElementPtr(pStatStreamRefNode);

      std::string statName
	= utl::dom::mustGetAttrString(pStatStreamRefElt,
				      eltName::statStreamRef_nameAttr);

      // For now, I'll just case this out on the statName.  These seem to
      // be living in the general pit of dumpables;
      //
      // ??????
      // It looks like there is nothing to prevent name collisions
      // among dumpables in the old code?
      // ??????
      //
      // An alternative might be to just create these.  Need
      // new memory management strategy for these.

      // This will lead us to the same exception if a bad name
      // (as per the schema) is used or if the dumpable is not
      // in the global pit.
      fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable = 0;
      if((statName == eltName::statStream_simTime)
	 || (statName == eltName::statStream_clockTime)
	 || (statName == eltName::statStream_speciesCount)
	 || (statName == eltName::statStream_reactionCount)
	 || (statName == eltName::statStream_reactionEventCount)
	 || (statName == eltName::statStream_volume))
	{
	  pDumpable = rMzrUnit.findDumpable(statName);
	}

      if(0 == pDumpable)
	throw unkStatStreamXcpt(statName,
				pStatStreamRefElt);

      return pDumpable;
    }
  };

  // This combinator class is necessitated by peculiar STL facts:
  // std::vector::push_back takes a reference argument, and this reference
  // argument type is used inappropriately by std::bind1st. (?)
  // Why do "adaptor" classes always seem to work so poorly for me?
  class addDumpableToTabDumpEvent :
    public std::unary_function<fnd::dumpable<fnd::basicDumpable::dumpArg>*, void>
  {
    tabDumpEvent* pEvent;
  public:
    addDumpableToTabDumpEvent(tabDumpEvent* pTabDumpEvent) :
      pEvent(pTabDumpEvent)
    {}

    void
    operator()(fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable) const
    {
      std::ostream& rOstream = pEvent->getOstream();
      
      fnd::dmpColumn<fnd::basicDumpable::dumpArg>* pColumn
	= new fnd::dmpColumn<fnd::basicDumpable::dumpArg>
	(pDumpable,
	 fnd::basicDumpable::dumpArg(rOstream));
      
      pEvent->push_back(pColumn);
    }
  };

  class parseDumpStream : public
  std::unary_function<xmlpp::Node*,
			  tabDumpEvent*>
  {
    mzrUnit& rMzrUnit;
    
  public:
    parseDumpStream(mzrUnit& refMzrUnit) :
      rMzrUnit(refMzrUnit)
    {}

    tabDumpEvent*
    operator()(xmlpp::Node* pDumpStreamNode) throw(std::exception)
    {
      xmlpp::Element* pDumpStreamElt
	= utl::dom::mustBeElementPtr(pDumpStreamNode);

      double dumpPeriod
	= utl::dom::mustGetAttrPosDouble
	(pDumpStreamElt,
	 eltName::dumpStream_dumpPeriodAttr);

      xmlpp::Element* pTargetFileElt
	= utl::dom::mustGetUniqueChild(pDumpStreamElt,
				       eltName::targetFile);
      std::string targetFileName
	= utl::dom::mustGetAttrString(pTargetFileElt,
				      eltName::targetFile_fileNameAttr);

      tabDumpEvent* pDumpEvent = new tabDumpEvent(dumpPeriod,
						  targetFileName);

      // Put the time in the first column.
      // It has to be looked up in the catalog of dumpables.
      std::string simTimeDumpableName("sim-time");
      fnd::dumpable<fnd::basicDumpable::dumpArg>* pSimTimeDumpable
	= rMzrUnit.mustFindDumpable(simTimeDumpableName,
				    pDumpStreamElt);
      std::ostream& rOstream = pDumpEvent->getOstream();

      fnd::dmpColumn<fnd::basicDumpable::dumpArg>* pColumn
	= new fnd::dmpColumn<fnd::basicDumpable::dumpArg>
	(pSimTimeDumpable,
	 fnd::basicDumpable::dumpArg(rOstream));
      
      pDumpEvent->push_back(pColumn);

      // Now we parse the different kinds of things that can be dumped.
      //
      // The speciesStreams are special in the role that they play
      // in dumping state.

      // Parse the species streams to get a vector of speciesDumpables.
      xmlpp::Node::NodeList speciesStreamRefNodes
	= pDumpStreamElt->get_children(eltName::speciesStreamRef);

      std::vector<fnd::dumpable<fnd::basicDumpable::dumpArg>*>
	speciesStreamDumpables(speciesStreamRefNodes.size());

      std::transform(speciesStreamRefNodes.begin(),
		     speciesStreamRefNodes.end(),
		     speciesStreamDumpables.begin(),
		     resolveSpeciesStreamRef(rMzrUnit));

      // Add these speciesDumpables to the tabDumpEvent.
      std::for_each(speciesStreamDumpables.begin(),
		    speciesStreamDumpables.end(),
		    addDumpableToTabDumpEvent(pDumpEvent));

      // Handle all the individual species.  These don't have
      // dumpables associated to them at this point.
      xmlpp::Node::NodeList speciesRefNodes
	= pDumpStreamElt->get_children(eltName::speciesRef);
      std::vector<singleSpeciesDumpable<mzrSpecies>*> singleSpeciesDumpables
	(speciesRefNodes.size());
      std::transform(speciesRefNodes.begin(),
		     speciesRefNodes.end(),
		     singleSpeciesDumpables.begin(),
		     resolveSpeciesRef(rMzrUnit));

      // Add these speciesDumpables to the tabDumpEvent.
      std::for_each(singleSpeciesDumpables.begin(),
		    singleSpeciesDumpables.end(),
		    addDumpableToTabDumpEvent(pDumpEvent));

      // Handle the stat streams.  These have dumpables in the global
      // pit, added in moleculizer::moleculizer().
      xmlpp::Node::NodeList statStreamRefNodes
	= pDumpStreamElt->get_children(eltName::statStreamRef);
      std::vector<fnd::dumpable<fnd::basicDumpable::dumpArg>*> statStreamDumpables
	(statStreamRefNodes.size());
      std::transform(statStreamRefNodes.begin(),
		     statStreamRefNodes.end(),
		     statStreamDumpables.begin(),
		     resolveStatStreamRef(rMzrUnit));

      // Add these dumpables (not speciesDumpables; statStreams don't do
      // anything special in a state dump.
      std::for_each(statStreamDumpables.begin(),
		    statStreamDumpables.end(),
		    addDumpableToTabDumpEvent(pDumpEvent));

      return pDumpEvent;
    }
  };

  class parseCreateEvent : public
  std::unary_function<xmlpp::Node*, void>
  {
    moleculizer& rMolzer;
    mzrUnit& rMzrUnit;
    
  public:
    parseCreateEvent(moleculizer& rMoleculizer,
		     mzrUnit& refMzrUnit) :
      rMolzer(rMoleculizer),
      rMzrUnit(refMzrUnit)
    {}
    
    void
    operator()(xmlpp::Node* pCreateEventNode) throw(std::exception)
    {
      xmlpp::Element* pCreateEventElt
	= utl::dom::mustBeElementPtr(pCreateEventNode);

      xmlpp::Element* pSpeciesRefElt
	= utl::dom::mustGetUniqueChild(pCreateEventElt,
				       eltName::speciesRef);
      std::string speciesName
	= utl::dom::mustGetAttrString(pSpeciesRefElt,
				      eltName::speciesRef_nameAttr);

      mzrSpecies* pSpecies
	= rMzrUnit.findSpecies(speciesName);
      if(0 == pSpecies) throw unkSpeciesXcpt(speciesName,
					     pSpeciesRefElt);

      // Get the number of molecules to create.
      xmlpp::Element* pPopulationElt
	= utl::dom::mustGetUniqueChild(pCreateEventElt,
				       eltName::population);
      int count
	= utl::dom::mustGetAttrPosInt(pPopulationElt,
				      eltName::population_countAttr);

      // Construct the event and install it in its deletion pit.
      createEvent* pCreateEvent = new createEvent(pSpecies,
						  count,
						  rMzrUnit);
      rMzrUnit.addUserEvent(pCreateEvent);

      // Get the time at which to schedule the event.
      // 
      // Note that we will probably want to be able to schedule events
      // at negative times.
      double time
	= utl::dom::mustGetAttrDouble(pCreateEventElt,
				      eltName::createEvent_timeAttr);

      // Schedule the event if it is supposed to occur in the future.
      // (which might not be the case if this simulation is being
      // restarted from a state dump.)
      double now = rMolzer.eventQ.getSimTime();
      if(now <= time)
	{
	  rMolzer.eventQ.scheduleEvent(pCreateEvent,
				       time);
	}
    }
  };


  class parseStopEvent : public
  std::unary_function<xmlpp::Node*, void>
  {
    moleculizer& rMolzer;
    mzrUnit& rMzrUnit;
    
  public:
    parseStopEvent(moleculizer& rMoleculizer,
		   mzrUnit& refMzrUnit) :
      rMolzer(rMoleculizer),
      rMzrUnit(refMzrUnit)
    {}
    
    void
    operator()(xmlpp::Node* pStopEventNode) const throw(std::exception)
    {
      xmlpp::Element* pStopEventElt
	= utl::dom::mustBeElementPtr(pStopEventNode);

      // Construct the event and install it in its deletion pit.
      stopEvent* pStopEvent = new stopEvent();
      rMzrUnit.addUserEvent(pStopEvent);

      // Note that we will in general want to be able to schedule
      // events at negative times (though probably not stop events.)
      double time
	= utl::dom::mustGetAttrDouble(pStopEventElt,
				      eltName::stopEvent_timeAttr);

      // Check that the event is supposed to happen in the future.
      double now = rMolzer.eventQ.getSimTime();
      if(time < now)
	throw stopEventInPastXcpt(now,
				  time);

      // Schedule the event
      rMolzer.eventQ.scheduleEvent(pStopEvent,
				    time);
    }
  };

  class parseDumpStateEvent :
    public std::unary_function<xmlpp::Node*, void>
  {
    moleculizer& rMolzer;
    mzrUnit& rMzrUnit;
    
  public:
    parseDumpStateEvent(moleculizer& rMoleculizer,
			mzrUnit& refMzrUnit) :
      rMolzer(rMoleculizer),
      rMzrUnit(refMzrUnit)
    {}
    
    void
    operator()(xmlpp::Node* pDumpStateEventNode) const throw(std::exception)
    {
      xmlpp::Element* pDumpStateEventElt
	= utl::dom::mustBeElementPtr(pDumpStateEventNode);

      // Get the target file for the state dump.
      xmlpp::Element* pTargetFileElt
	= utl::dom::mustGetUniqueChild(pDumpStateEventElt,
				       eltName::targetFile);

      std::string fileName
	= utl::dom::mustGetAttrString(pTargetFileElt,
				      eltName::targetFile_fileNameAttr);

      // Construct the dump event.
      dumpStateEvent* pDumpStateEvent
	= new dumpStateEvent(fileName);
      rMzrUnit.addUserEvent(pDumpStateEvent);

      // Parse the time at which the event should occur, and schedule the
      // event if it is in the future.  (When continuing a simulation
      // from a state dump, this gets rid of events in the moleculizer-input
      // that have already happened.
      double time
	= utl::dom::mustGetAttrDouble(pDumpStateEventElt,
				      eltName::dumpStateEvent_timeAttr);
      double now
	= rMolzer.eventQ.getSimTime();

      if(now <= time)
	{
	  rMolzer.eventQ.scheduleEvent(pDumpStateEvent,
				       time);
	}
    }
  };

  void
  mzrUnit::parseDomInput(xmlpp::Element* pRootElement,
			 xmlpp::Element* pModelElement,
			 xmlpp::Element* pStreamsElement,
			 xmlpp::Element* pEventsElement) throw(std::exception)
  {
    // Set the volume.  This is peculiar, and interesting, because
    // it's similar to setting the population of a species up front,
    // instead of running create events.  That's much more delicate
    // than this, and it's probably something that I won't let users
    // do.
    xmlpp::Element* pVolumeElt
      = utl::dom::mustGetUniqueChild(pModelElement,
				     eltName::volume);
    double volume
      = utl::dom::mustGetAttrDouble(pVolumeElt,
				    eltName::volume_litersAttr);
    // This is a bit funky: just discard the reactions that need updating,
    // since none of them are scheduled to begin with anyway?
    fnd::sensitivityList<mzrReaction> affectedReactions;
    getMolarFactor().updateVolume(volume,
				  affectedReactions);
  
    //////////////////////////////////////////////////////////////////

    // Get the explicit-reactions element.
    //
    // Eventually, there could be more than one kind of reaction,
    // so that this header element for all the different kinds isn't
    // useless.
    xmlpp::Element* pExplicitReactionsElt
      = utl::dom::mustGetUniqueChild(pModelElement,
				     eltName::explicitReactions);

    // Get the list of reaction elements.
    xmlpp::Node::NodeList reactionNodes
      = pExplicitReactionsElt->get_children(eltName::reaction);

    // Install each reaction in moleculizer's autoVector.
    std::for_each(reactionNodes.begin(),
		  reactionNodes.end(),
		  installReaction(*this));

    //////////////////////////////////////////////////////////////////

    // Get the dump-streams element under the streams element.
    //
    // At this point, there is only one kind of dump-stream, so the
    // dump-streams element looks redundant.
    xmlpp::Element* pDumpStreamsElt
      = utl::dom::mustGetUniqueChild(pStreamsElement,
				     eltName::dumpStreams);

    xmlpp::Node::NodeList dumpStreamNodes
      = pDumpStreamsElt->get_children(eltName::dumpStream);

    // Parse the dumpStream nodes to get the vector of tabDumpEvents,
    // are special among events in being needed when state is dumped.
    //
    // The header lines of dumpStreams are no longer printed as the
    // dumpStream input is "parsed."  Rather, they all print their
    // header lines in mzrUnit::prepareToRun.
    std::transform(dumpStreamNodes.begin(),
		   dumpStreamNodes.end(),
		   std::back_inserter(tabDumpEvents),
		   parseDumpStream(*this));

    // Insert the tabDumpEvents into the list of all events for management.
    std::copy(tabDumpEvents.begin(),
	      tabDumpEvents.end(),
	      std::back_inserter(userEvents));

    ///////////////////////////////////////////////////////////////////

    // Parse create events.
    xmlpp::Node::NodeList createEventNodes
      = pEventsElement->get_children(eltName::createEvent);
    std::for_each(createEventNodes.begin(),
		  createEventNodes.end(),
		  parseCreateEvent(rMolzer,
				   *this));

    // Parse stop-events.  I guess there will be only one, but let's
    // not assume so.
    xmlpp::Node::NodeList stopEventNodes
      = pEventsElement->get_children(eltName::stopEvent);

    // Issue a warning if there are no stop events.
    if(0 == stopEventNodes.size())
      noStopEventWarning(pEventsElement).warn();

    // Parse the stop events.
    std::for_each(stopEventNodes.begin(),
		  stopEventNodes.end(),
		  parseStopEvent(rMolzer,
				 *this));

    // Parse dump-state-events.
    xmlpp::Node::NodeList dumpStateEventNodes
      = pEventsElement->get_children(eltName::dumpStateEvent);
    std::for_each(dumpStateEventNodes.begin(),
		  dumpStateEventNodes.end(),
		  parseDumpStateEvent(rMolzer,
				      *this));
  }
}
