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

#include "cpt/cptEltName.hh"
#include "cpt/parseDumpStream.hh"
#include "cpt/singleGlobalSpeciesDumpable.hh"
#include "cpt/unkStatStreamXcpt.hh"
#include "cpt/unkGlobalDumpableXcpt.hh"

namespace cpt
{
  // For parsing a list of compartment-ref nodes, which gives a list of
  // compartments whose populations we want to dump, either totalled or
  // in separate columns.
  class parseDumpCompartment :
    public std::unary_function<xmlpp::Node*, void>
  {
    const compartmentGraph& rGraph;
    std::vector<int>& rCompartments;
    
  public:
    parseDumpCompartment(const compartmentGraph& rCompartmentGraph,
			 std::vector<int>& rDumpCompartments) :
      rGraph(rCompartmentGraph),
      rCompartments(rDumpCompartments)
    {}

    void
    operator()(xmlpp::Node* pCompartmentRefNode) const
    {
      xmlpp::Element* pCompartmentRefElt
	= utl::dom::mustBeElementPtr(pCompartmentRefNode);

      std::string compartmentName
	= utl::dom::mustGetAttrString(pCompartmentRefElt,
				      eltName::compartmentRef_nameAttr);

      int compartmentNdx
	= rGraph.mustFindCompartmentIndex(compartmentName,
					  pCompartmentRefElt);

      rCompartments.push_back(compartmentNdx);
    }
  };
  
  // Looks up dumpable for global species and sets up column(s) in dump
  // file for the compartment populations or the total population.
  class resolveSpeciesRef : public
  std::unary_function<xmlpp::Node*, void>
  {
    cptUnit& rCptUnit;
    tabDumpEvent* pEvent;
    
  public:
    resolveSpeciesRef(cptUnit& refCptUnit,
		      tabDumpEvent* pTabDumpEvent) :
      rCptUnit(refCptUnit),
      pEvent(pTabDumpEvent)
    {}

    void
    operator()(xmlpp::Node* pSpeciesRefNode) const
      throw(utl::xcpt)
    {
      xmlpp::Element* pSpeciesRefElt
	= utl::dom::mustBeElementPtr(pSpeciesRefNode);

      std::string speciesName
	= utl::dom::mustGetAttrString(pSpeciesRefElt,
				      eltName::speciesRef_nameAttr);

      // Look up the dumpable, which should have automatically
      // been created when the species was added to the cptUnit.
      fnd::basicDumpable*
	pBasicDumpable = rCptUnit.mustFindDumpable(speciesName,
						   pSpeciesRefElt);

      // This dynamic cast stinks yet again.  It's casting to the other main
      // case this time.
      fnd::dumpable<globalDumpArg>* pDumpable
	= dynamic_cast<fnd::dumpable<globalDumpArg>*>(pBasicDumpable);

      if(0 == pDumpable)
	throw unkGlobalDumpableXcpt(speciesName,
				    pSpeciesRefElt);

      // Determine whether or not compartment populations are to be summed.
      std::string totalPopsValue
	= utl::dom::mustGetAttrString(pSpeciesRefElt,
				      eltName::speciesRef_totalPopsAttr);
      bool sumCompartmentPops
	= (eltName::speciesRef_yesTotalPops == totalPopsValue);

      // If compartments are specified, parse them.  Otherwise, assume
      // we want to dump (or total) all compartments.
      xmlpp::Element* pDumpCompartmentsElt
	= utl::dom::getOptionalChild(pSpeciesRefElt,
				     eltName::dumpCompartments);
      std::vector<int> dumpCompartmentIndices;
      if(pDumpCompartmentsElt)
	{
	  xmlpp::Node::NodeList compartmentRefNodes
	    = pDumpCompartmentsElt->get_children(eltName::compartmentRef);
	  std::for_each(compartmentRefNodes.begin(),
			compartmentRefNodes.end(),
			parseDumpCompartment(rCptUnit.getCompartmentGraph(),
					     dumpCompartmentIndices));
	}
      else
	{
	  int compartmentCount
	    = rCptUnit.getCompartmentGraph().compartments.size();
	  
	  for(int dumpCompartmentNdx = 0;
	      dumpCompartmentNdx < compartmentCount;
	      ++dumpCompartmentNdx)
	    {
	      dumpCompartmentIndices.push_back(dumpCompartmentNdx);
	    }
	}
      
      // Push the dumpable, together with the right kind of dumpArg,
      // onto the tabDumpEvent's columns.  For now, the column
      // should give the total of all the compartment populations.
      //
      // The tabDumpEvent does memory management of these.
      fnd::dmpColumn<globalDumpArg>* pColumn
	= new fnd::dmpColumn<globalDumpArg>
	(pDumpable,
	 globalDumpArg(pEvent->getOstream(),
		       dumpCompartmentIndices,
		       sumCompartmentPops));
      pEvent->push_back(pColumn);
    }
  };

  // Looks up dumpables connected with omniplexes, etc, and adds them
  // to a tabDumpEvent being parsed.
  class resolveSpeciesStreamRef : public
  std::unary_function<xmlpp::Node*, void>
  {
    cptUnit& rCptUnit;
    tabDumpEvent* pEvent;

  public:
    resolveSpeciesStreamRef(cptUnit& refCptUnit,
			    tabDumpEvent* pTabDumpEvent) :
      rCptUnit(refCptUnit),
      pEvent(pTabDumpEvent)
    {}

    void
    operator()(xmlpp::Node* pSpeciesStreamRefNode) const
      throw(std::exception)
    {
      xmlpp::Element* pSpeciesStreamRefElt
	= utl::dom::mustBeElementPtr(pSpeciesStreamRefNode);

      std::string dumpableName
	= utl::dom::mustGetAttrString(pSpeciesStreamRefElt,
				      eltName::speciesStreamRef_nameAttr);

      // Look up the speciesDumpable.
      fnd::basicDumpable* pBasicDumpable
	= rCptUnit.findDumpable(dumpableName);
      
      // This dynamic cast stinks yet again.  It's casting to the other main
      // case this time.
      fnd::dumpable<globalDumpArg>* pDumpable
	= dynamic_cast<fnd::dumpable<globalDumpArg>*>(pBasicDumpable);

      if(0 == pDumpable)
	throw unkGlobalDumpableXcpt(dumpableName,
				    pSpeciesStreamRefElt);

      // Determine whether or not compartment populations are to be summed.
      std::string totalPopsValue
	= utl::dom::mustGetAttrString(pSpeciesStreamRefElt,
				      eltName::speciesStreamRef_totalPopsAttr);
      bool sumCompartmentPops
	= (eltName::speciesStreamRef_yesTotalPops == totalPopsValue);

      // If compartments are specified, parse them.  Otherwise, assume
      // we want to dump (or total) all compartments.
      xmlpp::Element* pDumpCompartmentsElt
	= utl::dom::getOptionalChild(pSpeciesStreamRefElt,
				     eltName::dumpCompartments);
      std::vector<int> dumpCompartmentIndices;
      if(pDumpCompartmentsElt)
	{
	  xmlpp::Node::NodeList compartmentRefNodes
	    = pDumpCompartmentsElt->get_children(eltName::compartmentRef);
	  std::for_each(compartmentRefNodes.begin(),
			compartmentRefNodes.end(),
			parseDumpCompartment(rCptUnit.getCompartmentGraph(),
					     dumpCompartmentIndices));
	}
      else
	{
	  int compartmentCount
	    = rCptUnit.getCompartmentGraph().compartments.size();
	  
	  for(int dumpCompartmentNdx = 0;
	      dumpCompartmentNdx < compartmentCount;
	      ++dumpCompartmentNdx)
	    {
	      dumpCompartmentIndices.push_back(dumpCompartmentNdx);
	    }
	}
      
      // Push the dumpable, together with the right kind of dumpArg,
      // onto the tabDumpEvent's columns.  For now, the column
      // should give the total of all the compartment populations.
      //
      // The tabDumpEvent does memory management of these.
      fnd::dmpColumn<globalDumpArg>* pColumn
	= new fnd::dmpColumn<globalDumpArg>
	(pDumpable,
	 globalDumpArg(pEvent->getOstream(),
		       dumpCompartmentIndices,
		       sumCompartmentPops));
      pEvent->push_back(pColumn);
    }
  };

  class resolveStatStreamRef : public
  std::unary_function<xmlpp::Node*, void>
  {
    cptUnit& rCptUnit;
    tabDumpEvent* pEvent;

  public:
    resolveStatStreamRef(cptUnit& refCptUnit,
			 tabDumpEvent* pTabDumpEvent) :
      rCptUnit(refCptUnit),
      pEvent(pTabDumpEvent)
    {}

    void
    operator()(xmlpp::Node* pStatStreamRefNode) const
      throw(utl::xcpt)
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
      fnd::basicDumpable* pBasicDumpable = 0;
      if((statName == eltName::statStream_simTime)
	 || (statName == eltName::statStream_clockTime)
	 || (statName == eltName::statStream_speciesCount)
	 || (statName == eltName::statStream_reactionCount)
	 || (statName == eltName::statStream_reactionEventCount)
	 || (statName == eltName::statStream_volume))
	{
	  pBasicDumpable = rCptUnit.findDumpable(statName);
	}
      
      // Make sure the dumpable is of the right type.  This stinks,
      // somewhat.
      fnd::dumpable<fnd::basicDumpable::dumpArg>* pDumpable
	= dynamic_cast<fnd::dumpable<fnd::basicDumpable::dumpArg>*>(pBasicDumpable);

      // Throw unknown dumpable exception (normally thrown by
      // cptUnit::mustFindDumpable) if the name was bad or if the
      // dumpable was not of the right type.  This looks like a source
      // of confusion.
      if(0 == pDumpable)
	throw unkStatStreamXcpt(statName,
				pStatStreamRefElt);

      // The columns are memory managed by the dumpable that they go into.
      fnd::dmpColumn<fnd::basicDumpable::dumpArg>* pColumn
	= new fnd::dmpColumn<fnd::basicDumpable::dumpArg>
	(pDumpable,
	 fnd::basicDumpable::dumpArg(pEvent->getOstream()));

      pEvent->push_back(pColumn);
    }
  };

  void
  parseDumpStream::
  operator()(xmlpp::Node* pDumpStreamNode)
    throw(utl::xcpt)
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

    // Create the tabDumpEvent, and register it in the cptUnit for memory
    // mangagement.
    //
    // This also adds this event to the tabDumpEvents vector, which is used
    // when state is dumped to emit the tags of the species whose populations
    // are totalled in the various columns.
    //
    // tabDumpEvents are initialized and scheduled for the first time
    // in cptUnit::prepareToRun.
    tabDumpEvent* pDumpEvent = new tabDumpEvent(dumpPeriod,
						targetFileName);
    rCptUnit.addTabDumpEvent(pDumpEvent);

    // Put the time in the first column.
    // It has to be looked up in the catalog of dumpables.
    std::string simTimeDumpableName("sim-time");
    fnd::basicDumpable* pBasicDumpable
      = rCptUnit.mustFindDumpable(simTimeDumpableName,
				  pDumpStreamElt);
    // This cast stinks, again!
    fnd::dumpable<fnd::basicDumpable::dumpArg>* pSimTimeDumpable
      = dynamic_cast<fnd::dumpable<fnd::basicDumpable::dumpArg>*>
      (pBasicDumpable);

    std::ostream& rOstream = pDumpEvent->getOstream();

    // Memory management is done by the tabDumpEvent.
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

    // Look up the associated "global" dumpables and push them
    // onto the dumpEvent's vector with the right kind of dumpArg.
    std::for_each(speciesStreamRefNodes.begin(),
		  speciesStreamRefNodes.end(),
		  resolveSpeciesStreamRef(rCptUnit,
					  pDumpEvent));

    // Handle all the individual species.
    xmlpp::Node::NodeList speciesRefNodes
      = pDumpStreamElt->get_children(eltName::speciesRef);

    std::for_each(speciesRefNodes.begin(),
		  speciesRefNodes.end(),
		  resolveSpeciesRef(rCptUnit,
				    pDumpEvent));

    // Handle the stat streams.  These have dumpables in the global
    // pit, added in moleculizer::moleculizer().
    xmlpp::Node::NodeList statStreamRefNodes
      = pDumpStreamElt->get_children(eltName::statStreamRef);

    std::for_each(statStreamRefNodes.begin(),
		  statStreamRefNodes.end(),
		  resolveStatStreamRef(rCptUnit,
				       pDumpEvent));
  }
}
