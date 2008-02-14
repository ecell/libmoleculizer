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

#include "cpt/cptApp.hh"
#include "cpt/respondReaction.hh"
#include "cpt/createEvent.hh"
#include "cpt/parseGlobalSpecies.hh"
#include "clx/cptPlexFamily.hh"
#include "clx/parsePlex.hh"
#include "clx/parseOmniPlex.hh"

namespace clx
{
  namespace
  {
    // Class for accumulating list of all plex nodes for omniplexes.
    class addPathNodesToList :
      public std::unary_function<std::string, void>
    {
      xmlpp::Node::NodeList& rNodes;
      xmlpp::Element* pRootElt;
    public:
      addPathNodesToList(xmlpp::Node::NodeList& rNodeList,
			 xmlpp::Element* pRootElement) :
	rNodes(rNodeList),
	pRootElt(pRootElement)
      {}
  
      void
      operator()(const std::string& rXpath) const
	throw()
      {
	// Get the list satisfying this xPath query.
	xmlpp::NodeSet xpathHits
	  = pRootElt->find(rXpath);

	// Insert these at the end of the list of all hits.
	rNodes.insert(rNodes.end(),
		      xpathHits.begin(),
		      xpathHits.end());
      }
    };
  }

  void
  clxUnit::parseDomInput(xmlpp::Element* pRootElt,
			  xmlpp::Element* pModelElt,
			  xmlpp::Element* pStreamsElt,
			  xmlpp::Element* pEventsElt)
    throw(utl::xcpt)
  {
    // First we pick out a number of elements and lists of elements
    // for tranversal.

    // Species streams.
    xmlpp::Element* pSpeciesStreamsElt
      = utl::dom::mustGetUniqueChild(pStreamsElt,
				     cpt::eltName::speciesStreams);

    xmlpp::Node::NodeList omniSpeciesStreamNodes
      = pSpeciesStreamsElt
      ->get_children(eltName::omniSpeciesStream);

    xmlpp::Node::NodeList plexSpeciesStreamNodes
      = pSpeciesStreamsElt
      ->get_children(eltName::plexSpeciesStream);

    // Allosteric omnis.
    xmlpp::Element* pAlloOmnisElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     eltName::allostericOmnis);
    xmlpp::Node::NodeList alloOmniNodes
      = pAlloOmnisElt->get_children(eltName::allostericOmni);

    // Allosteric plexes.
    xmlpp::Element* pAlloPlexesElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     eltName::allostericPlexes);
    xmlpp::Node::NodeList alloPlexNodes
      = pAlloPlexesElt->get_children(eltName::allostericPlex);

    // Explicit plexSpecies.
    xmlpp::Element* pExplicitSpeciesElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     cpt::eltName::explicitSpecies);
    xmlpp::Node::NodeList plexSpeciesNodes
      = pExplicitSpeciesElt->get_children(eltName::plexSpecies);

    // Use Xpaths to omniplex nodes, registered by other modules,
    // to find all the omniplexes in the file.
    //
    // We have to have all the omniplexes in place before we recognize any
    // complexes in the conventional way, thereby creating plexFamilies.
    xmlpp::Node::NodeList omniPlexNodes;
    std::for_each(omniXpaths.begin(),
		  omniXpaths.end(),
		  addPathNodesToList(omniPlexNodes,
				     pRootElt));

    // "Unify" omniplex families (recognize, but without the usual
    // initializations) and put them on the clxUnit's list of omniplexes.
    // After this is done, plexes and omniplexes can be recognized in the usual
    // way.
    std::for_each(omniPlexNodes.begin(),
		  omniPlexNodes.end(),
		  parseOmniPlex(rCptUnit,
				rCmlUnit,
				*this));

    // Since the omniplex families have been "unified" in, they won't
    // undergo normal initialization when recognized.  Hence, we
    // run through them all and connect them to their features.
    //
    // After this point, we should be ready to ready to recognize complexes
    // in the usual way, creating plexFamilies.
    std::for_each(omniPlexFamilies.begin(),
		  omniPlexFamilies.end(),
		  std::mem_fun(&cptPlexFamily::connectToFeatures));

    // Parse allosteric omnis.
    //
    // This must be completed before any species of complexes are generated.
    std::for_each(alloOmniNodes.begin(),
		  alloOmniNodes.end(),
		  parseAllostericOmni(rCmlUnit,
				      *this));

    // Parse allosteric plexes.
    //
    // This must be done before any species of complexes are generated.
    std::for_each(alloPlexNodes.begin(),
		  alloPlexNodes.end(),
		  parseAllostericPlex(rCmlUnit,
				      *this,
				      rCptUnit));
		  

    // Attach dumpables to families of complexes.  This must be done before
    // any species of complexes are generated.

    // Parse query-based dumpables for omniplexes.
    std::for_each(omniSpeciesStreamNodes.begin(),
		  omniSpeciesStreamNodes.end(),
		  parseOmniSpeciesStream(rCptUnit,
					 rCmlUnit,
					 *this));

    // Parse query-based dumpables for plexes.
    std::for_each(plexSpeciesStreamNodes.begin(),
		  plexSpeciesStreamNodes.end(),
		  parsePlexSpeciesStream(rCptUnit,
					 rCmlUnit,
					 *this));

    // Parse explicit plexSpecies, generating species of complexes, but not
    // populating them.  Since this doesn't create the initial population,
    // notification isn't an issue.
    std::for_each(plexSpeciesNodes.begin(),
		  plexSpeciesNodes.end(),
		  parseExplicitPlexSpecies(rCptUnit,
					   rCmlUnit,
					   *this));
  }

  namespace
  {
    // Class for creating the initial populations of explicit plex species.
    //
    // This is the first time that notification happens, at least for
    // plexSpecies.
    class createInitialPop : public
    std::unary_function<xmlpp::Node*, void>
    {
      cpt::cptApp& rCptApp;
      cpt::cptUnit& rCptUnit;
    public:
      createInitialPop(cpt::cptApp& rCptApp,
		       cpt::cptUnit& rCpt) :
	rCptApp(rCptApp),
	rCptUnit(rCpt)
      {}
    
      void
      operator()(xmlpp::Node* pPlexSpeciesNode) const
	throw(utl::xcpt)
      {
	xmlpp::Element* pPlexSpeciesElt
	  = utl::dom::mustBeElementPtr(pPlexSpeciesNode);

	std::string speciesName
	  = utl::dom::mustGetAttrString(pPlexSpeciesElt,
					eltName::plexSpecies_nameAttr);

	// Parse the initial populations of the species.
	std::vector<int> compartmentPops
	  = cpt::parseCompartmentPops(pPlexSpeciesNode,
				      rCptUnit.getCompartmentGraph());
	cpt::globalSpecies* pSpecies
	  = rCptUnit.mustFindSpecies(speciesName,
				     pPlexSpeciesElt);

	// Rather than scheduling an event, we effectively do a create event.
	cpt::createEvent creator(pSpecies,
				 compartmentPops);
	creator.happen(rCptApp);
      }
    };
  }

  void
  clxUnit::
  prepareToRun(xmlpp::Element* pRootElt,
	       xmlpp::Element* pModelElt,
	       xmlpp::Element* pStreamsElt,
	       xmlpp::Element* pEventsElt)
    throw(utl::xcpt)
  {
    // Create the initial population of all explicit plex species.
    xmlpp::Element* pExplicitSpeciesElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     cpt::eltName::explicitSpecies);
    xmlpp::Node::NodeList plexSpeciesNodes
      = pExplicitSpeciesElt->get_children(eltName::plexSpecies);

    std::for_each(plexSpeciesNodes.begin(),
		  plexSpeciesNodes.end(),
		  createInitialPop(rCptApp,
				   rCptUnit));
  }

  namespace
  {
    class accumulateSpecies : public
    std::unary_function<std::multimap<int, cptPlexFamily*>::value_type, void>
    {
      std::vector<cptPlexSpecies*>& rAllSpecies;
    public:
      accumulateSpecies(std::vector<cptPlexSpecies*>& rAllSpeciesVector) :
	rAllSpecies(rAllSpeciesVector)
      {}

      void
      operator()(const argument_type& rEntry) const
      {
	cptPlexFamily* pPlexFamily = rEntry.second;
	pPlexFamily->accumulateSpecies(rAllSpecies);
      }
    };

    class zeroUpdateSpecies : public
    std::unary_function<cptPlexSpecies*, void>
    {
      fnd::sensitivityList<cpt::cptReaction>& rAffected;
      std::vector<int> compartmentDeltas;
      int depth;
    
    public:
      zeroUpdateSpecies(fnd::sensitivityList<cpt::cptReaction>& rAffectedReactions,
			int compartmentCount,
			int notifyDepth) :
	rAffected(rAffectedReactions),
	compartmentDeltas(compartmentCount, 0)
      {}

      void
      operator()(cptPlexSpecies* pPlexSpecies) const
      {
	pPlexSpecies->update(compartmentDeltas,
			     rAffected,
			     depth);
      }
    };
  }

  void
  clxUnit::prepareToDump(xmlpp::Element* pRootElt,
			  xmlpp::Element* pModelElt,
			  xmlpp::Element* pStreamsElt,
			  xmlpp::Element* pEventsElt,
			  xmlpp::Element* pTaggedSpeciesElement)
    throw(utl::xcpt)
  {
    // Parse the tagged-plex-species just about like explicit plex-species,
    // just ignoring the tag.  The population will also end up being ingnored.
    // (For explicit plex-species, they are traversed again, in order to
    // initialize the populations.)
    xmlpp::Node::NodeList taggedPlexSpeciesNodes
      = pTaggedSpeciesElement->get_children(eltName::taggedPlexSpecies);

    // Going along with the now-well-established principle of stupid names.
    // The analogous routine for reading explicit plex-species is
    // "processPlexSpecies".
    std::vector<cptPlexSpecies*> updatedSpecies;
    std::for_each(taggedPlexSpeciesNodes.begin(),
		  taggedPlexSpeciesNodes.end(),
		  parseTaggedPlexSpecies(rCptUnit,
					 rCmlUnit,
					 *this,
					 updatedSpecies));

    // Run through all the plexSpecies, and update each by 0.
    // This should cause all the reactions for the species to come
    // into existence.  This may be a few more reactions than
    // were dumped, since some of the dumped species (with population 0)
    // may not ever have been updated.


    // Update each plexSpecies by zero, accumulating the
    // set of affected reactions.
    int compartmentCount
      = (int) rCptUnit.getCompartmentGraph().compartments.size();
    fnd::sensitivityList<cpt::cptReaction> affectedReactions;
    std::for_each(updatedSpecies.begin(),
		  updatedSpecies.end(),
		  zeroUpdateSpecies(affectedReactions,
				    compartmentCount,
				    rCptUnit.getGenerateDepth()));

    // Checking on just how bad the inflation of species and reactions is.
    std::vector<cptPlexSpecies*> afterPlexSpecies;
    std::for_each(recognize.plexHasher.begin(),
		  recognize.plexHasher.end(),
		  accumulateSpecies(afterPlexSpecies));

    // For consistency's sake, rescheduling the affected reactions.
    // Since the update was by 0, this isn't really necessary.
    std::for_each(affectedReactions.begin(),
		  affectedReactions.end(),
		  cpt::respondReaction(rCptApp.getPropensities()));

    // This is unit's default implementation of prepareToDump.
    prepareToRun(pRootElt,
		 pModelElt,
		 pStreamsElt,
		 pEventsElt);
  }

  class parseTaggedPlexSpeciesNpop :
    public std::unary_function<xmlpp::Node*, void>
  {
    cpt::cptUnit& rCptUnit;
    cml::cmlUnit& rCmlUnit;
    clxUnit& rClxUnit;
    fnd::sensitivityList<cpt::cptReaction>& rAffected;
    int depth;
    
  public:
    parseTaggedPlexSpeciesNpop
    (cpt::cptUnit& refCptUnit,
     cml::cmlUnit& refCmlUnit,
     clxUnit& refClxUnit,
     fnd::sensitivityList<cpt::cptReaction>& rAffectedReactions,
     int notifyDepth) :
      rCptUnit(refCptUnit),
      rCmlUnit(refCmlUnit),
      rClxUnit(refClxUnit),
      rAffected(rAffectedReactions),
      depth(notifyDepth)
    {}
    
    void
    operator()(xmlpp::Node* pTaggedPlexSpeciesNode) const
      throw(utl::xcpt)
    {
      xmlpp::Element* pTaggedPlexSpeciesElt
	= utl::dom::mustBeElementPtr(pTaggedPlexSpeciesNode);

      // In this case, we expect at most one species to appear
      // in updatedSpecies; in the usual application of
      // parseTaggedPlexSpecies, many species are parsed,
      // and the ones that had been updated at the time of dump
      // appear in updatedSpecies.
      std::vector<cptPlexSpecies*> updatedSpecies;
      parseTaggedPlexSpecies tpsParser(rCptUnit,
				       rCmlUnit,
				       rClxUnit,
				       updatedSpecies);
      cptPlexSpecies* pSpecies
	= tpsParser(pTaggedPlexSpeciesElt);
      

      // Parse the initial populations of the species.
      std::vector<int> compartmentPops
	= cpt::parseCompartmentPops(pTaggedPlexSpeciesElt,
				    rCptUnit.getCompartmentGraph());

      // Update the species with its dumped population.  This obviates
      // the "zeroUpdate" pass that has to be made in prepareToDump,
      // as well as the usual "prepareToRun" which updates the
      // explicit plex species with the populations provided in
      // cptApp-input.
      xmlpp::Node::NodeList updatedNodes
	= pTaggedPlexSpeciesElt->get_children(eltName::updated);

      // Update the species if called for.
      //
      // The plexSpecies we just parsed above is the only one that could
      // appear in updatedSpecies, due to this special application of
      // parseTaggedPlexSpecies. In the usual application of
      // parseTaggedPlexSpecies, many species are parsed, and the ones that
      // had been updated at the time of dump appear in updatedSpecies.
      if(0 < updatedSpecies.size())
	{
	  pSpecies->update(compartmentPops,
			   rAffected,
			   depth);
	}
    }
  };

  void
  clxUnit::prepareToContinue(xmlpp::Element* pRootElt,
			      xmlpp::Element* pModelElt,
			      xmlpp::Element* pStreamsElt,
			      xmlpp::Element* pEventsElt,
			      std::map<std::string, std::string>& rTagToName,
			      xmlpp::Element* pTaggedSpeciesElement)
    throw(utl::xcpt)
  {
    // Parse the tagged-plex-species.
    xmlpp::Node::NodeList taggedPlexSpeciesNodes
      = pTaggedSpeciesElement->get_children(eltName::taggedPlexSpecies);

    // Going along with the now-well-established principle of stupid names.
    // The analogous routine for reading explicit plex-species is
    // "processPlexSpecies".
    //
    // In this version the populations from the dump are used as the
    // initial populations of the species.  All the plex species that
    // exist should be updated here, so that the "zeroUpdate" pass
    // is not necessary.
    fnd::sensitivityList<cpt::cptReaction> affectedReactions;
    std::for_each(taggedPlexSpeciesNodes.begin(),
		  taggedPlexSpeciesNodes.end(),
		  parseTaggedPlexSpeciesNpop(rCptUnit,
					     rCmlUnit,
					     *this,
					     affectedReactions,
					     rCptUnit.getGenerateDepth()));

    // Reschedule the affected reactions.  The current simulation time
    // should have been set to the dump time by the cptUnit.
    std::for_each(affectedReactions.begin(),
		  affectedReactions.end(),
		  cpt::respondReaction(rCptApp.getPropensities()));

    // In this version, do NOT run prepareToRun here, as that
    // sets the populations of the explicit species and reschedules
    // reactions.  May want rearchitecture here.
  }
}
