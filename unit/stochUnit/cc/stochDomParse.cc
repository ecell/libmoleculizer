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

#include "mzr/moleculizer.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/mzrUnitParse.hh"
#include "mol/molEltName.hh"
#include "stoch/stochSpecies.hh"
#include "stoch/stochUnit.hh"
#include "stoch/stochEltName.hh"

namespace stoch
{
  class unnamedStochSpeciesTagXcpt : public mzr::mzrXcpt
  {
    std::string
    mkMsg(xmlpp::Node* pOffendingNode)
    {
      std::ostringstream msgStream;
      msgStream << domUtils::domXcptMsg(pOffendingNode)
		<< "Stoch species tag not named in explicit-species-tags.";
      return msgStream.str();
    }
  public:
    unnamedStochSpeciesTagXcpt(xmlpp::Node* pOffendingNode) :
      mzr::mzrXcpt(mkMsg(pOffendingNode))
    {}
  };
	
  class parseStochSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    stochUnit& rStochUnit;
    
  public:
    parseStochSpecies(stochUnit& refStochUnit) :
      rStochUnit(refStochUnit)
    {}
    
    void
    operator()(xmlpp::Node* pStochSpeciesNode) const throw(std::exception)
    {
      xmlpp::Element* pStochSpeciesElt
	= domUtils::mustBeElementPtr(pStochSpeciesNode);

      // Get the name of the stochastirator species.
      std::string speciesName
	= domUtils::mustGetAttrString
	(pStochSpeciesElt,
	 eltName::stochSpecies_nameAttr);

      // Get the molecular weight of the stochastirator species.
      xmlpp::Element* pWeightElt
	= domUtils::mustGetUniqueChild(pStochSpeciesElt,
				       bnd::eltName::weight);
      double molWeight
	= domUtils::mustGetAttrPosDouble
	(pWeightElt,
	 bnd::eltName::weight_daltonsAttr);

      // Construct the stochastirator species.
      stochSpecies* pSpecies = new stochSpecies(speciesName,
						molWeight);

      // Add to the catalog of stoch species.
      rStochUnit.addStochSpecies(speciesName,
				 pSpecies);
    }
  };

  void
  stochUnit::parseDomInput(xmlpp::Element* pRootElement,
			   xmlpp::Element* pModelElement,
			   xmlpp::Element* pStreamsElement,
			   xmlpp::Element* pEventsElement) throw(std::exception)
  {
    // Parse all the stoch species as named species.
    xmlpp::Element* pExplicitSpeciesElt
      = domUtils::mustGetUniqueChild(pModelElement,
				     mzr::eltName::explicitSpecies);

    xmlpp::Node::NodeList stochSpeciesNodes
      = pExplicitSpeciesElt->get_children(eltName::stochSpecies);

    std::for_each(stochSpeciesNodes.begin(),
		  stochSpeciesNodes.end(),
		  parseStochSpecies(*this));
  }

  class createInitialPop :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::moleculizer& rMolzer;
    mzr::mzrUnit& rMzrUnit;
  public:
    createInitialPop(mzr::moleculizer& rMoleculizer,
		     mzr::mzrUnit& rMzr) :
      rMolzer(rMoleculizer),
      rMzrUnit(rMzr)
    {
    }

    void
    operator()(xmlpp::Node* pStochSpeciesNode) const
      throw(std::exception)
    {
      xmlpp::Element* pStochSpeciesElt
	= domUtils::mustBeElementPtr(pStochSpeciesNode);

      std::string speciesName
	= domUtils::mustGetAttrString(pStochSpeciesElt,
				      eltName::stochSpecies_nameAttr);

      xmlpp::Element* pPopulationElt
	= domUtils::mustGetUniqueChild(pStochSpeciesElt,
				       mzr::eltName::population);
      int pop
	= domUtils::mustGetAttrInt(pPopulationElt,
				   mzr::eltName::population_countAttr);
      mzr::species* pSpecies
	= rMzrUnit.findSpecies(speciesName);

      if(0 == pSpecies)
	throw mzr::unknownSpeciesXcpt(pStochSpeciesElt,
				      speciesName);

      mzr::createEvent creator(pSpecies,
			       pop,
			       rMzrUnit);
      creator.doEvent(rMolzer,
		      rMzrUnit);
    }
  };

  void
  stochUnit::prepareToRun(xmlpp::Element* pRootElt,
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
      = domUtils::mustGetUniqueChild(pModelElt,
				     mzr::eltName::explicitSpecies);

    xmlpp::Node::NodeList stochSpeciesNodes
      = pExplicitSpeciesElt->get_children(eltName::stochSpecies);

    std::for_each(stochSpeciesNodes.begin(),
		  stochSpeciesNodes.end(),
		  createInitialPop(rMolzer,
				   rMzrUnit));
  }

  class processTaggedStochSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    std::set<mzr::reaction*>& rAffectedReactions;
    std::map<std::string, std::string>& rTagToName;
    mzr::mzrUnit& rMzr;
    
  public:
    processTaggedStochSpecies
    (std::set<mzr::reaction*>& refAffectedReactions,
     std::map<std::string, std::string>& rTagToNameMap,
     mzr::mzrUnit& rMzrUnit) :
      rAffectedReactions(refAffectedReactions),
      rTagToName(rTagToNameMap),
      rMzr(rMzrUnit)
    {}

    void
    operator()(xmlpp::Node* pTaggedStochSpeciesNode) const
      throw(std::exception)
    {
      xmlpp::Element* pTaggedStochSpeciesElt
	= domUtils::mustBeElementPtr(pTaggedStochSpeciesNode);

      std::string tag
	= domUtils::mustGetAttrString(pTaggedStochSpeciesElt,
				      eltName::taggedStochSpecies_tagAttr);

      // Look up the species name in the tag to name map.
      std::map<std::string, std::string>::iterator iEntry
	= rTagToName.find(tag);

      if(rTagToName.end() == iEntry)
	throw unnamedStochSpeciesTagXcpt(pTaggedStochSpeciesElt);

      std::string speciesName = iEntry->second;

      // Look up the species in mzrUnit's catalog.
      mzr::species* pSpecies = rMzr.findSpecies(speciesName);
      if(0 == pSpecies)
	throw mzr::unknownSpeciesXcpt(pTaggedStochSpeciesElt,
				      speciesName);


      // Get the population of the stoch species when dumped.
      xmlpp::Element* pPopulationElt
	= domUtils::mustGetUniqueChild(pTaggedStochSpeciesElt,
				       eltName::population);

      int pop
	= domUtils::mustGetAttrInt(pPopulationElt,
				   eltName::population_countAttr);

      // Update the species.
      pSpecies->update(rAffectedReactions,
		       pop);
    }
  };

  void
  stochUnit::prepareToContinue(xmlpp::Element* pRootElt,
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
    std::set<mzr::reaction*> affectedReactions;
    std::for_each(taggedStochSpeciesNodes.begin(),
		  taggedStochSpeciesNodes.end(),
		  processTaggedStochSpecies(affectedReactions,
					    rTagToName,
					    rMzrUnit));

    // Reschedule the affected reactions.
    std::for_each(affectedReactions.begin(),
		  affectedReactions.end(),
		  mzr::scheduleReaction(rMolzer.eventQ,
					rMzrUnit));

    // In this version, do NOT run prepareToRun, as that sets the
    // populations of the stochSpecies to that given in the explicit
    // stoch-species elements of moleculizer-input.
  }
}
