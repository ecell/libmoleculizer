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
#include "cpt/dupSpeciesNameXcpt.hh"
#include "cpt/unitsMgr.hh"
#include "cpt/parseGlobalSpecies.hh"
#include "cpt/respondReaction.hh"
#include "cst/unkStochSpeciesXcpt.hh"
#include "cst/badStochSpeciesTagXcpt.hh"
#include "cst/cstEltName.hh"
#include "cst/cstUnit.hh"


namespace cst
{
  cstUnit::
  cstUnit(cpt::cptApp& rCptApp,
	  cpt::cptUnit& refCptUnit) :
    cpt::unit("cst",
	      rCptApp),
    rCptUnit(refCptUnit)
  {
    // Unit isn't responsible for any model elements or
    // reaction generators.

    // Add responsibility for stoch species.
    inputCap.addExplicitSpeciesContentName(eltName::stochSpecies);

    // Not responsible for any species streams or events.
  }

  bool
  cstUnit::
  addStochSpecies(cStochSpecies* pStochSpecies)
  {
    // Add the species to the overall catalog of named species,
    // creating a dumpable for it at the same time.
    bool result = rCptUnit.addNamedSpecies(pStochSpecies->getName(),
					   pStochSpecies);

    // Add to memory-management in this unit.
    return result && stochSpeciesByName.addEntry(pStochSpecies->getName(),
						 pStochSpecies);
  }

  void
  cstUnit::
  mustAddStochSpecies(cStochSpecies* pStochSpecies,
		      xmlpp::Node* pRequestingNode)
    throw(utl::xcpt)
  {
    if(! addStochSpecies(pStochSpecies))
      throw cpt::dupSpeciesNameXcpt(pStochSpecies->getName(),
				    pRequestingNode);
  }
  
  cStochSpecies*
  cstUnit::
  findStochSpecies(const std::string& rSpeciesName) const
  {
    return stochSpeciesByName.findEntry(rSpeciesName);
  }

  cStochSpecies*
  cstUnit::
  mustFindStochSpecies(const std::string& rSpeciesName,
		       xmlpp::Node* pRequestingNode) const
    throw(utl::xcpt)
  {
    cStochSpecies* pSpecies = findStochSpecies(rSpeciesName);

    if(0 == pSpecies) throw unkStochSpeciesXcpt(rSpeciesName,
						pRequestingNode);
    return pSpecies;
  }

  void
  cstUnit::
  addNoSubstrateArrow(cpt::cptReaction* pReaction)
  {
    noSubstrateArrows.push_back(pReaction);
  }

  class insertStochSpeciesElt :
    public std::unary_function<utl::autoCatalog<cStochSpecies>::value_type, void>
  {
    xmlpp::Element* pTaggedSpeciesElt;
    
  public:
    insertStochSpeciesElt(xmlpp::Element* pTaggedSpeciesElement) :
      pTaggedSpeciesElt(pTaggedSpeciesElement)
    {}

    void
    operator()(const argument_type& rCatalogEntry) const
      throw(std::exception)
    {
      const cStochSpecies* pStochSpecies = rCatalogEntry.second;

      pStochSpecies->insertElt(pTaggedSpeciesElt);
    }
  };

  void
  cstUnit::
  insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
  {
    xmlpp::Element* pModelElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     cpt::eltName::model);

    // Ensure that the tagged-species element is there.
    xmlpp::Element* pTaggedSpeciesElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     cpt::eltName::taggedSpecies);

    // Insert stoch-species nodes.
    std::for_each(stochSpeciesByName.begin(),
		  stochSpeciesByName.end(),
		  insertStochSpeciesElt(pTaggedSpeciesElt));

    // Leaving out the "no-substrate reactions" for now.  Don't know
    // why I put them in this unit to begin with anyway...
  }

  class parseCompartmentPop :
    public std::unary_function<xmlpp::Node*, std::pair<int, int> >
  {
    const cpt::compartmentGraph& rGraph;
    
  public:
    parseCompartmentPop(const cpt::compartmentGraph& rCompartmentGraph) :
      rGraph(rCompartmentGraph)
    {}
    
    std::pair<int, int>
    operator()(xmlpp::Node* pCompartmentPopNode) const
    {
      xmlpp::Element* pCompartmentPopElt
	= utl::dom::mustBeElementPtr(pCompartmentPopNode);

      std::string compartmentName
	= utl::dom::mustGetAttrString(pCompartmentPopElt,
				      eltName::compartmentPop_compartmentNameAttr);
      int compartmentNdx
	= rGraph.mustFindCompartmentIndex(compartmentName,
					  pCompartmentPopElt);

      int pop
	= utl::dom::mustGetAttrNNInt(pCompartmentPopElt,
				     eltName::compartmentPop_populationAttr);

      return std::pair<int, int>(compartmentNdx, pop);
    }
  };

  class createInitialPop :
    public std::unary_function<xmlpp::Node*, void>
  {
    cpt::cptApp& rApp;
    const cstUnit& rUnit;
    const cpt::compartmentGraph& rGraph;
    
  public:
    createInitialPop(cpt::cptApp& rCptApp,
		     cstUnit& rCstUnit) :
      rApp(rCptApp),
      rUnit(rCstUnit),
      rGraph(rApp.getUnits()->getCptUnit().getCompartmentGraph())
    {}

    void
    operator()(xmlpp::Node* pStochSpeciesNode) const
      throw(std::exception)
    {
      xmlpp::Element* pStochSpeciesElt
	= utl::dom::mustBeElementPtr(pStochSpeciesNode);

      // Locate the species by its name.
      std::string speciesName
	= utl::dom::mustGetAttrString(pStochSpeciesElt,
				      eltName::stochSpecies_nameAttr);
      cStochSpecies* pSpecies
	= rUnit.mustFindStochSpecies(speciesName,
				     pStochSpeciesElt);

      // Parse the initial populations of the species.
      std::vector<int> compartmentPops
	= cpt::parseCompartmentPops(pStochSpeciesNode,
				    rGraph);

      // Update the species, collecting the reactions whose propensities
      // need to be changed because of the population changes.
      fnd::sensitivityList<cpt::cptReaction> affectedReactions;
      pSpecies->update(compartmentPops,
		       affectedReactions,
		       rApp.getCptUnit().getGenerateDepth());

      // Update the propensities of the affected reactions.
      std::for_each(affectedReactions.begin(),
		    affectedReactions.end(),
		    cpt::respondReaction(rApp.getPropensities()));
    }
  };

  void
  cstUnit::prepareToRun(xmlpp::Element* pRootElt,
			  xmlpp::Element* pModelElt,
			  xmlpp::Element* pStreamsElt,
			  xmlpp::Element* pEventsElt)
    throw(std::exception)
  {
    // Create the initial population of all explicit stochSpecies.
    //
    // stoch::prepareToRun is also stochUnit's prepareToDump; this
    // works because stochSpecies don't participate in automatic
    // reaction generation (so all stochSpecies are explicit as are
    // all the reactions in which they participate.)
    xmlpp::Element* pExplicitSpeciesElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     cpt::eltName::explicitSpecies);

    xmlpp::Node::NodeList stochSpeciesNodes
      = pExplicitSpeciesElt->get_children(eltName::stochSpecies);

    std::for_each(stochSpeciesNodes.begin(),
		  stochSpeciesNodes.end(),
		  createInitialPop(rCptApp,
				   *this));
  }

  class processTaggedStochSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    std::map<std::string, std::string>& rTagToName;
    const cstUnit& rUnit;
    cpt::cptApp& rApp;
    const cpt::compartmentGraph& rGraph;
    
  public:
    processTaggedStochSpecies(std::map<std::string, std::string>& rTagToNameMap,
			      cstUnit& rCstUnit,
			      cpt::cptApp& rCptApp) :
      rTagToName(rTagToNameMap),
      rUnit(rCstUnit),
      rApp(rCptApp),
      rGraph(rCptApp.getUnits()->getCptUnit().getCompartmentGraph())
    {}

    void
    operator()(xmlpp::Node* pTaggedStochSpeciesNode) const
      throw(std::exception)
    {
      xmlpp::Element* pTaggedStochSpeciesElt
	= utl::dom::mustBeElementPtr(pTaggedStochSpeciesNode);

      std::string tag
	= utl::dom::mustGetAttrString(pTaggedStochSpeciesElt,
				      eltName::taggedStochSpecies_tagAttr);

      // Look up the species name in the tag to name map.
      std::map<std::string, std::string>::iterator iEntry
	= rTagToName.find(tag);

      if(rTagToName.end() == iEntry)
	throw badStochSpeciesTagXcpt(pTaggedStochSpeciesElt);

      std::string speciesName = iEntry->second;

      // Look up the species in the catalog.
      cStochSpecies* pSpecies
	= rUnit.mustFindStochSpecies(speciesName,
				     pTaggedStochSpeciesElt);

      // Parse the compartment populations.
      std::vector<int> compartmentPops
	= cpt::parseCompartmentPops(pTaggedStochSpeciesElt,
				    rGraph);

      // Update the species, collecting the reactions whose propensities
      // need to be changed because of the population changes.
      fnd::sensitivityList<cpt::cptReaction> affectedReactions;
      pSpecies->update(compartmentPops,
		       affectedReactions,
		       rApp.getCptUnit().getGenerateDepth());

      // Update the propensities of the affected reactions.
      std::for_each(affectedReactions.begin(),
		    affectedReactions.end(),
		    cpt::respondReaction(rApp.getPropensities()));
    }
  };

  void
  cstUnit::
  prepareToContinue(xmlpp::Element* pRootElt,
		    xmlpp::Element* pModelElt,
		    xmlpp::Element* pStreamsElt,
		    xmlpp::Element* pEventsElt,
		    std::map<std::string, std::string>& rTagToName,
		    xmlpp::Element* pTaggedSpeciesElement)
    throw(std::exception)
  {
    // Parse the tagged-stoch-species.
    xmlpp::Node::NodeList taggedStochSpeciesNodes
      = pTaggedSpeciesElement->get_children(eltName::taggedStochSpecies);

    // This will also update each stoch species with the population
    // given in the state dump.
    std::for_each(taggedStochSpeciesNodes.begin(),
		  taggedStochSpeciesNodes.end(),
		  processTaggedStochSpecies(rTagToName,
					    *this,
					    rCptApp));

    // In this version, do NOT run prepareToRun, as that sets the
    // populations of the stochSpecies to that given in the explicit
    // stoch-species elements of moleculizer-input.
  }
}
