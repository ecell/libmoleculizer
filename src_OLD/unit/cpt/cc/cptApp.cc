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

#include <limits>
#include <fstream>
#include <iostream>
#include "utl/forBoth.hh"
#include "utl/badNNIntArgXcpt.hh"
#include "utl/unkArgXcpt.hh"
#include "cpt/respondReaction.hh"
#include "cpt/parseCompartmentGraph.hh"
#include "cpt/parseGlobalSpecies.hh"
#include "cpt/parseGlobalReaction.hh"
#include "cpt/unitsMgr.hh"
#include "cpt/cptUnit.hh"
#include "cpt/inputCapTest.hh"
#include "cpt/strayModelContentXcpt.hh"
#include "cpt/strayReactionGenXcpt.hh"
#include "cpt/strayExplicitSpeciesContentXcpt.hh"
#include "cpt/straySpeciesStreamContentXcpt.hh"
#include "cpt/strayEventContentXcpt.hh"

namespace cpt
{
  void
  cptApp::
  constructorPrelude(void)
  {
    // Add up the input capabilities of the units.
    pUnits->unionInputCaps(inputCap);

    // Formerly initialized random seed here.
  }

  // Used in contructorCore below.
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



  void
  cptApp::
  constructorCore(xmlpp::Element* pRootElement,
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
      throw strayModelContentXcpt(*iUnhandledModelContent);

    // Get the reaction-gens node, which contains kinds of reaction generators
    // introduced by units.
    xmlpp::Element* pReactionGensElt
      = utl::dom::mustGetUniqueChild(pModelElement,
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
      throw strayReactionGenXcpt(*iUnhandledReactionGenContent);

    // Get the explicit species model node, which can contain kinds of explicit
    // species introduced by units.
    xmlpp::Element* pExplicitSpeciesElt
      = utl::dom::mustGetUniqueChild(pModelElement,
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
	strayExplicitSpeciesContentXcpt(*iUnhandledExplicitSpeciesContent);

    // Get the speciesStreams node, which can contain kinds of species streams
    // introduced by units.
    xmlpp::Element* pSpeciesStreamsElt
      = utl::dom::mustGetUniqueChild(pStreamsElement,
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
      throw straySpeciesStreamContentXcpt(*iUnhandledSpeciesStreamsContent);
  
    // Verify that every element in the events section is handled by some
    // unit or another.
    xmlpp::Node::NodeList eventsContentNodes
      = pEventsElement->get_children();
    xmlpp::Node::NodeList::iterator iUnhandledEventsContent
      = std::find_if(eventsContentNodes.begin(),
		     eventsContentNodes.end(),
		     eventNodeNotInCap(inputCap));
    if(eventsContentNodes.end() != iUnhandledEventsContent)
      throw strayEventContentXcpt(*iUnhandledEventsContent);

    // Parse the compartment graph.
    //
    // This was in cptUnit, but it needs to be done before any unit
    // can parse, e.g. globalSpecies.
    xmlpp::Element* pCompartmentGraphElt
      = utl::dom::mustGetUniqueChild(pModelElement,
				     eltName::compartmentGraph);
    parseCompartmentGraph graphParser(*this);
    graphParser(pCompartmentGraphElt);

    // Have each unit do its parsing thing.
    std::for_each(pUnits->begin(),
		  pUnits->end(),
		  unitParseDomInput(pRootElement,
				    pModelElement,
				    pStreamsElement,
				    pEventsElement));
  }

  class prepareUnitToRun :
    public std::unary_function<unit*, void>
  {
    xmlpp::Element* pRootElt;
    xmlpp::Element* pModelElt;
    xmlpp::Element* pStreamsElt;
    xmlpp::Element* pEventsElt;
  public:
    prepareUnitToRun(xmlpp::Element* pRootElement,
		     xmlpp::Element* pModelElement,
		     xmlpp::Element* pStreamsElement,
		     xmlpp::Element* pEventsElement) :
      pRootElt(pRootElement),
      pModelElt(pModelElement),
      pStreamsElt(pStreamsElement),
      pEventsElt(pEventsElement)
    {
    }

    void
    operator()(unit* pUnit) const
      throw(std::exception)
    {
      pUnit->prepareToRun(pRootElt,
			  pModelElt,
			  pStreamsElt,
			  pEventsElt);
    }
  };

  cptApp::
  cptApp(int argCount,
	 char** argVector,
	 xmlpp::Document* pDoc)
    throw(std::exception) :
    pUnits(new unitsMgr(*this)),
    rCptUnit(pUnits->getCptUnit()),
    now(0.0),
    propensities(rCptUnit.getRng())
  {
    // Do the "input capabilities" thing.
    constructorPrelude();

    // Get command-line arguments.
    rCptUnit.parseCommandLine(argCount,
			      argVector);
    
    xmlpp::Element* pRootElement
      = pDoc->get_root_node();

    // Parse run parameters.
    xmlpp::Element* pRunParamsElt
      = utl::dom::mustGetUniqueChild(pRootElement,
				     eltName::runParams);
    cycleTime
      = utl::dom::mustGetAttrPosDouble(pRunParamsElt,
				       eltName::runParams_cycleTimeAttr);
    reactionsPerCycle
      = utl::dom::mustGetAttrPosDouble(pRunParamsElt,
				       eltName::runParams_reactionsPerCycleAttr);
    epsilon
      = utl::dom::mustGetAttrPosDouble(pRunParamsElt,
				       eltName::runParams_epsilonAttr);

    // Get the unique model element.
    xmlpp::Element* pModelElement
      = utl::dom::mustGetUniqueChild(pRootElement,
				     eltName::model);

    // Get the unique streams element.
    xmlpp::Element* pStreamsElement
      = utl::dom::mustGetUniqueChild(pRootElement,
				     eltName::streams);

    // Get the unique events element.
    xmlpp::Element* pEventsElement
      = utl::dom::mustGetUniqueChild(pRootElement,
				     eltName::events);

    // Extract model info.
    constructorCore(pRootElement,
		    pModelElement,
		    pStreamsElement,
		    pEventsElement);
  
    // Have each unit do its prepareToRun thing.
    std::for_each(pUnits->begin(),
		  pUnits->end(),
		  prepareUnitToRun(pRootElement,
				   pModelElement,
				   pStreamsElement,
				   pEventsElement));
  }

  cptApp::
  ~cptApp(void)
  {
    delete pUnits;
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

  // State dump, invoked in a dump-state event.
  xmlpp::Document*
  cptApp::makeDomOutput(void)
    throw(std::exception)
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
    std::for_each(getUnits()->begin(),
		  getUnits()->end(),
		  unitInsertStateElements(pRootElt));

    return pDoc;
  }

  int
  cptApp::
  run(void)
    throw(std::exception)
  {
    // Construct the multinomial sampler.
    //
    // This is done here because we are now sure that the compartmentGraph
    // has its final number of compartments.
    utl::gsl::multinomialSampler sampler
      (rCptUnit.getCompartmentGraph().compartments.size(),
       rCptUnit.getRng());

    int cycleCount = 0;
    int eventCount = 0;
    simulationCycle(sampler,
		    cycleCount,
		    eventCount);

    std::cerr << "Run consisted of "
	      << cycleCount
	      << " cycles, and entailed "
	      << eventCount
	      << " reaction events."
	      << std::endl;

    return 0;
  }

  // Supporting determination of the time step.
  class checkDiffusionMatrixNorm :
    public std::unary_function<globalSpecies*, void>
  {
    double& rMaxNorm;

  public:
    checkDiffusionMatrixNorm(double& rMaxDiffusionMatrixNorm) :
      rMaxNorm(rMaxDiffusionMatrixNorm)
    {}

    void
    operator()(const globalSpecies* pSpecies) const
    {
      double speciesDiffusionNorm
	= pSpecies->getDiffusionMatrix().boxNorm();

      if(rMaxNorm < speciesDiffusionNorm) rMaxNorm = speciesDiffusionNorm;
    }
  };

  // Supporting the diffusion phase of the simulation cycle.
  class diffuseSpecies :
    public std::unary_function<globalSpecies*, void>
  {
    utl::gsl::multinomialSampler& rSampler;
    double cycleTime;
    double epsilon;
    fnd::sensitivityList<cptReaction>& rAffected;
    int depth;
    
  public:
    diffuseSpecies(utl::gsl::multinomialSampler& rMultinomialSampler,
		   double theCycleTime,
		   double theEpsilon,
		   fnd::sensitivityList<cptReaction>& rAffectedReactions,
		   int generateDepth) :
      rSampler(rMultinomialSampler),
      cycleTime(theCycleTime),
      epsilon(theEpsilon),
      rAffected(rAffectedReactions),
      depth(generateDepth)
    {}

    void
    operator()(globalSpecies* pGlobalSpecies) const
    {
      pGlobalSpecies->doDiffusion(rSampler,
				  cycleTime,
				  epsilon,
				  rAffected,
				  depth);
    }

  };

  // This routine exists mainly so that I can escape the simulation cycle
  // neatly (by returning from this routine) without throwing/catching my way
  // out.
  void
  cptApp::
  simulationCycle(utl::gsl::multinomialSampler& rSampler,
		  int& rCycleCount,
		  int& rEventCount)
  {
    // This involves sweeping over all global species in a couple of places.
    const std::vector<globalSpecies*>& rAllGlobalSpecies
      = rCptUnit.getGlobalSpecies();
    
    // Construct Poisson sampler for determining the (maximum) number of
    // reactions in a cycle, given the total propensity of all the reactions
    // and the time step.
    utl::gsl::poissonSampler poisson(rCptUnit.getRng());

    // Get the maximum of the norms of the diffusion matrices for all the
    // global species, for aid in controlling the length of the time step.
    double maxDiffusionNorm = 0;
    std::for_each(rAllGlobalSpecies.begin(),
		  rAllGlobalSpecies.end(),
		  checkDiffusionMatrixNorm(maxDiffusionNorm));

    // Cook up limit on time step from factors not changing during simulation.
    double maxStep = cycleTime;
        
    while(true)
      {
	// Calculate the time step.

	// The time step should be no longer than the maximum step size.
	double step = maxStep;

	// The time step should be no longer than the time to the next user
	// event.
	double timeToNextEvent = getQueue().getNextEventTime() - getSimTime();
	if(timeToNextEvent < step) step = timeToNextEvent;

	// The time step should be short enough that we expect no more
	// reactions in it than given by the reactionsPerCycle global.
	double totalPropensity = getPropensities().getTotalPropensity();
	if(0.0 < totalPropensity)
	  {
	    double propensityStep
	      = reactionsPerCycle / totalPropensity;

	    if(propensityStep < step) step = propensityStep;
	  }

	// Sample Poisson with rate step * totalPropensity to get the number
	// of reactions for this cycle.
	//
	// Note that, much of the time, this cancels out the above calcuation,
	// and effectively samples Poisson with rate reactionsPerCycle.
	double poissonRate = step * totalPropensity;
	int cycleReactionCount = 0;
	if (0.0 < poissonRate)
	  {
	    cycleReactionCount = poisson.sample(poissonRate);
	  }

	// Sample propensity distro, do reactions, until we've done the number
	// of reactions calculated above, or until there are no reactions that
	// can be done.
	while(0 < cycleReactionCount)
	  {
	    cptReaction* pReaction = getPropensities().sample()->second;

	    pReaction->happen(*this);
	    ++rEventCount;

	    if(getPropensities().getTotalPropensity() <= 0.0)
	      {
		cycleReactionCount = 0;
	      }
	    else
	      {
		--cycleReactionCount;
	      }
	  }

	// Do diffusion.  Diffusion can cause new global species to be
	// generated and added to the allGlobalSpecies vector, so this
	// traversal can't be done with an iterator.
	fnd::sensitivityList<cptReaction> affectedReactions;
	int globalSpeciesNdx = 0;
	diffuseSpecies speciesDiffuser(rSampler,
				       step,
				       getEpsilon(),
				       affectedReactions,
				       -1);
	//				       rCptUnit.getGenerateDepth());
	while(globalSpeciesNdx < (int) rAllGlobalSpecies.size())
	  {
	    speciesDiffuser(rAllGlobalSpecies[globalSpeciesNdx]);

	    ++globalSpeciesNdx;
	  }
	
	// Adjust compartment reaction propensities, total propensity
	// after diffusion.
	std::for_each(affectedReactions.begin(),
		      affectedReactions.end(),
		      respondReaction(getPropensities()));

	// Do "pending" user events; i.e. those whose scheduled times
	// are before the end of this cycle.
	double endOfCycle = getSimTime() + step;
	double nextEventTime = theQueue.getNextEventTime();
	// Comparison with <= is necessary, since the next event
	// may happen precisely at the end of the cycle.
	while(nextEventTime <= endOfCycle)
	  {
	    // Get the next event, erasing its entry from the queue.
	    cptEvent* pEvent = theQueue.getNextEvent();

	    // Do the event at the time when it was scheduled to occur,
	    // possibly ending the simulation.
	    setSimTime(nextEventTime);
	    if(fnd::stop == pEvent->happen(*this)) return;

	    nextEventTime = theQueue.getNextEventTime();
	  }

	// Advance time to the end of the cycle.
	setSimTime(endOfCycle);
	++rCycleCount;
      }
  }
}
