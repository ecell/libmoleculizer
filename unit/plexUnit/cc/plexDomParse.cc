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
#include "mzr/mzrUnitParse.hh"
#include "mol/siteShape.hh"
#include "mol/molState.hh"
#include "mol/modMol.hh"
#include "mol/modQuery.hh"
#include "mol/molDomParse.hh"
#include "plex/plexEltName.hh"
#include "plex/prm.hh"
#include "plex/plexFamily.hh"
#include "plex/plexUnit.hh"
#include "plex/plexDomParse.hh"
#include "plex/parseOmniPlex.hh"

namespace plx
{

  // Adds a binding to a parser plex (using the parser plex's
  // instance-name ==> instance-index map.)
  class addBindingToPlex :
    public std::unary_function<xmlpp::Node*, void>
  {
    parserPlex& rPlex;
    // To avoid making more than one binding involving a given binding site.
    // This replaces a very inefficient way of doing this in the old scripting
    // approach.
    std::set<plexSiteSpec>& rBound;
  public:
    addBindingToPlex(parserPlex& rParserPlex,
		     std::set<plexSiteSpec>& rBoundSites) :
      rPlex(rParserPlex),
      rBound(rBoundSites)
    {}

    void
    operator()(xmlpp::Node* pBindingNode) const throw(mzr::mzrXcpt)

  };

  // This determines the plexFamily of an XML plex in a minimal way.  The
  // recognized plexFamily is not initialized in the usual way, and the
  // mapping from mol instance names to (paradigm) mol indexes is not
  // constructed.  This function is intended for initial processing of
  // omniplexes.
  plexFamily*
  unifyPlexElt(xmlpp::Element* pPlexElt,
	       bnd::molUnit& rMolUnit,
	       plexUnit& rPlexUnit)
  {
    parserPlex parsedPlex;
  
    // Process the mol instances.
    xmlpp::Node::NodeList molInstances
      = pPlexElt->get_children(eltName::molInstance);
    std::for_each(molInstances.begin(),
		  molInstances.end(),
		  addMolInstance(rMolUnit,
				 parsedPlex));

    // Process the bindings.
    //
    // This is to make sure that no site is in two bindings.
    std::set<plexSiteSpec> boundSites;
    xmlpp::Node::NodeList bindings
      = pPlexElt->get_children(eltName::binding);
    std::for_each(bindings.begin(),
		  bindings.end(),
		  addBindingToPlex(parsedPlex,
				   boundSites));

    // Recognize the plex, but without doing any of the usual
    // initializations of the resulting plexFamily.
    plexFamily* returnValue = 0;
    rPlexUnit.recognize.unify(parsedPlex,
			      returnValue);

    return returnValue;
  }

  // Used to "adapt" the instance name to instance index map of a
  // parserPlex to the indexing of the plexFamily's paradigm plex,
  // using the recognition isomorphism.
  class remapMolNdx :
    public std::unary_function<const std::pair<const std::string, int>&,
    std::pair<std::string, int> >
  {
    const std::vector<int>& rMolMap;
  public:
    remapMolNdx(const plexIsoPair& rIsoPair) :
      rMolMap(rIsoPair.forward.molMap)
    {}

    std::pair<std::string, int>
    operator()(const std::pair<const std::string, int>& rSourceEntry) const
    {
      return std::make_pair(rSourceEntry.first,
			    rMolMap[rSourceEntry.second]);
    }
  };

  // This determines the plexFamily of an XML plex in the usual way,
  // initializing the plexFamily if it is new, and generating a parserPlex
  // whose instance name to mol index mapping is "adapted" to the plexFamily's
  // paradigm.
  //
  // This parserPlex can then be used to deal with name-based definitions of
  // queries, allostery specifications, etc.
  plexFamily*
  recognizePlexElt(xmlpp::Element* pPlexElt,
		   parserPlex& rParserPlex,
		   bnd::molUnit& rMolUnit,
		   plexUnit& rPlexUnit)
  {
    parserPlex parsedPlex;
  
    // Process the mol instances.
    xmlpp::Node::NodeList molInstances
      = pPlexElt->get_children(eltName::molInstance);
    std::for_each(molInstances.begin(),
		  molInstances.end(),
		  addMolInstance(rMolUnit,
				 parsedPlex));

    // Process the bindings.
    //
    // This is to make sure that no site is in two bindings.
    std::set<plexSiteSpec> boundSites;
    xmlpp::Node::NodeList bindings
      = pPlexElt->get_children(eltName::binding);
    std::for_each(bindings.begin(),
		  bindings.end(),
		  addBindingToPlex(parsedPlex,
				   boundSites));

    // Recognize the plex, but without doing any of the usual
    // initializations of the resulting plexFamily.
    plexIsoPair recognitionIsos;
    plexFamily* pFamily = rPlexUnit.recognize(parsedPlex,
					      recognitionIsos);

    // Reindex the instance-name to instance-index map.
    std::map<std::string, int> nameToMolNdx;
    std::transform(parsedPlex.nameToMolNdx.begin(),
		   parsedPlex.nameToMolNdx.end(),
		   std::inserter(nameToMolNdx,
				 nameToMolNdx.begin()),
		   remapMolNdx(recognitionIsos));

    rParserPlex = parserPlex(nameToMolNdx,
			     pFamily->getParadigm());

    return pFamily;
  }

  std::pair<plexSiteSpec, bnd::siteParam>
  makeAlloSiteMapEntry::
  operator()(xmlpp::Node* pInstanceRefNode) const
    throw(mzr::mzrXcpt)
  {
    // Check that the node is an element node, possibly unnecessarily.
    xmlpp::Element* pInstanceRefElt
      = domUtils::mustBeElementPtr(pInstanceRefNode);
    
    // Get the mol instance name.
    std::string molInstanceName
      = domUtils::mustGetAttrString
      (pInstanceRefElt,
       eltName::molInstanceRef_nameAttr);

    // Get the corresponding mol index, part of the plexSiteSpec, in order to
    // make an entry in the allostery map of the plex definition.
    int molNdx = rPlex.getMolNdxByName(molInstanceName);
    if(0 > molNdx) throw unknownMolInstanceXcpt(pInstanceRefElt,
						molInstanceName);

    // Get the corresponding mol, in order to look up binding site by
    // name.
    bnd::mol* pMol = rPlex.mols[molNdx];

    // Get the child element that gives the site name.
    xmlpp::Element* pBindingSiteRefElt
      = domUtils::mustGetUniqueChild(pInstanceRefElt,
				     bnd::eltName::bindingSiteRef);

    // Get the name of the binding site.
    std::string bindingSiteName
      = domUtils::mustGetAttrString
      (pBindingSiteRefElt,
       bnd::eltName::bindingSiteRef_nameAttr);

    // Get the index of the binding site, part of the plexSiteSpec.
    int bindingSiteNdx = -1;
    if(! pMol->findSite(bindingSiteName,
			bindingSiteNdx))
      throw bnd::unknownSiteXcpt(pBindingSiteRefElt,
				 bindingSiteName);

    // Get the actual binding site, to look up the shape.
    bnd::bindingSite& rBindingSite = pMol->getSiteByNdx(bindingSiteNdx);
    

    // Get the child element that gives the shape of the binding
    // site when in this subcomplex.
    xmlpp::Element* pSiteShapeRefElt
      = domUtils::mustGetUniqueChild(pBindingSiteRefElt,
				     bnd::eltName::siteShapeRef);

    // Get the name of the site shape.
    std::string siteShapeName
      = domUtils::mustGetAttrString
      (pSiteShapeRefElt,
       bnd::eltName::siteShapeRef_nameAttr);

    // Look up the site shape.
    const bnd::siteShape* pSiteShape = rBindingSite.getParam(siteShapeName);
    if(0 == pSiteShape)
      throw bnd::unknownSiteShapeXcpt(pSiteShapeRefElt,
				      rBindingSite,
				      pMol,
				      siteShapeName);

    // Return entry for allostery map.
    return std::make_pair(plexSiteSpec(molNdx,
				       bindingSiteNdx),
			  pSiteShape);
  }

  // Processes a single modification map entry into a query for a plex dumpable
  // or omniplex dumpable.
  class addModClauseToQuery :
    public std::unary_function<xmlpp::Node*, void>
  {
    bnd::molUnit& rMolUnit;
    plexUnit& rPlexUnit;
    andPlexQueries* pQuery;
    bnd::modMol* pMol;
    plexMolSpec instanceNdx;
  public:
    addModClauseToQuery(bnd::molUnit& refMolUnit,
			plexUnit& refPlexUnit,
			andPlexQueries* pAndPlexQueries,
			bnd::modMol* pModMol,
			plexMolSpec molInstanceNdx) :
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      pQuery(pAndPlexQueries),
      pMol(pModMol),
      instanceNdx(molInstanceNdx)
    {}

    void
    operator()(xmlpp::Node* pModSiteRefNode) const throw(std::exception)
    {
      // Do the element cast, probably unnecessarily dynamic.
      xmlpp::Element* pModSiteRefElt
	= domUtils::mustBeElementPtr(pModSiteRefNode);
    
      // Get the modification site name.
      std::string modSiteName
	= domUtils::mustGetAttrString(pModSiteRefElt,
				      bnd::eltName::modSiteRef_nameAttr);

      // Look up the modification site on the mol, getting the
      // index of the modification site.
      int modSiteNdx = pMol->getModSiteNdx(modSiteName);
      if(0 > modSiteNdx)
	throw bnd::unknownModSiteXcpt(pModSiteRefElt,
				      modSiteName,
				      pMol);

      // Get the modification name.
      xmlpp::Element* pModRefElt
	= domUtils::mustGetUniqueChild(pModSiteRefNode,
				       bnd::eltName::modRef);
      std::string modName
	= domUtils::mustGetAttrString(pModRefElt,
				      bnd::eltName::modRef_nameAttr);

      // Look up the modification.
      const bnd::modification* pMod
	= rMolUnit.getMod(modName);
      if(0 == pMod)
	throw bnd::unknownModXcpt(pModRefElt,
				  modName);

      // Construct the query and register for memory management
      // by the plexUnit.
      bnd::modMixinStateQuery* pModQuery
	= new bnd::modMixinStateQuery(instanceNdx,
				 modSiteNdx,
				 pMod);
      rPlexUnit.addPlexQuery(pModQuery);

      // Add the query to the conjunction of queries to be applied to
      // species of complexes appearing in this stream.
      pQuery->addQuery(pModQuery);
    }
  };

  // Processes modification map entries into mol state queries, which
  // are added to a query for dumping.
  void
  addInstanceClausesToQuery::
  operator()(xmlpp::Node* pModMolInstanceRefNode) const
    throw(std::exception)
  {
    // Make sure the instance-state is an element, probably unnecessarily.
    xmlpp::Element* pModMolInstanceRefElt
      = dynamic_cast<xmlpp::Element*>(pModMolInstanceRefNode);
    if(0 == pModMolInstanceRefElt)
      throw domUtils::badElementCastXcpt(pModMolInstanceRefNode);

    std::string instanceName
      = domUtils::mustGetAttrString
      (pModMolInstanceRefElt,
       eltName::modMolInstanceRef_nameAttr);

    // Get the instance index.  Each filter clause needs a plexMolSpec,
    // which is an instance index together with the index of the binding
    // on the instance mol.
    int instanceNdx = rPlex.getMolNdxByName(instanceName);
    if(0 > instanceNdx)
      throw unknownMolInstanceXcpt(pModMolInstanceRefElt,
				   instanceName);


    bnd::mol* pMol = rPlex.mols[instanceNdx];
    bnd::modMol* pModMol = dynamic_cast<bnd::modMol*>(pMol);
    if(0 == pModMol)
      throw bnd::badModMolCastXcpt(pModMolInstanceRefElt,
				   pMol);


    // Process the modification map, which gives the state
    // of this instance mol.  This is a map from binding sites to
    // shapes.
    xmlpp::Element* pModMapElt
      = domUtils::mustGetUniqueChild(pModMolInstanceRefElt,
				     bnd::eltName::modMap);
    xmlpp::Node::NodeList modSiteRefs
      = pModMapElt->get_children(bnd::eltName::modSiteRef);
    std::for_each(modSiteRefs.begin(),
		  modSiteRefs.end(),
		  addModClauseToQuery(rMolUnit,
				      rPlexUnit,
				      pQuery,
				      pModMol,
				      instanceNdx));
  }

  // Installs allosteric modifications into plexFamily.  This can be applied
  // to an allosteric-plex node or an allosteric-omni node, which have the
  // same content (defined in the schema by "allo-plex-content").
  class doAlloMods : public
  std::unary_function<xmlpp::Node*, void>
  {
    bnd::molUnit& rMolUnit;
    plexUnit& rPlexUnit;
    
  public:
    doAlloMods(bnd::molUnit& refMolUnit,
	       plexUnit& refPlexUnit) :
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}
    
    void
    operator()(xmlpp::Node* pAlloPlexNode) const throw(mzr::mzrXcpt)
    {
      // Process the complex structure, and recognize it (in the usual,
      // family-initializing way.)
      xmlpp::Element* pPlexElt
	= domUtils::mustGetUniqueChild(pAlloPlexNode,
				       eltName::plex);
      parserPlex parsedPlex;
      plexFamily* pFamily = recognizePlexElt(pPlexElt,
					     parsedPlex,
					     rMolUnit,
					     rPlexUnit);

      // Process the query specifications: a list of mol instances
      // with a "state filter" for each.  Now we have made this element
      // optional, so that allosteric-plex and allosteric-omni can use
      // it without their having to be changed.
      xmlpp::Node::NodeList instanceStatesNodes
	= pAlloPlexNode->get_children(eltName::instanceStates);

      // Construct the query and register it for memory management
      // by the plexUnit.
      andPlexQueries* pPlexQueries = new andPlexQueries();
      rPlexUnit.addPlexQuery(pPlexQueries);

      // Process queries from the file.
      if(instanceStatesNodes.size() != 0)
	{
	  xmlpp::Element* pInstanceStatesElt
	    = domUtils::mustBeElementPtr(instanceStatesNodes.front());
	  
	  // Eventually we will have to parse state filters for different
	  // kinds of mols, but for now, just mod-mols
	  xmlpp::Node::NodeList modMolInstanceRefNodes
	    = pInstanceStatesElt->get_children
	    (eltName::modMolInstanceRef);

	  // Make queries out of the mod-mol instance specifications;
	  // all the queries are AND'ed together by the query
	  // pointed to by pPlexQueries.
	  std::for_each(modMolInstanceRefNodes.begin(),
			modMolInstanceRefNodes.end(),
			addInstanceClausesToQuery(pPlexQueries,
						  parsedPlex,
						  rMolUnit,
						  rPlexUnit));
	}


      // Process the allosteric sites.
      //
      // Is the "allosteric-sites" element dispensable?
      xmlpp::Element* pAlloSitesElt
	= domUtils::mustGetUniqueChild(pAlloPlexNode,
				       eltName::allostericSites);

      // Process allosteric sites.
      siteToShapeMap alloSiteMap;
      xmlpp::Node::NodeList instanceRefNodes
	= pAlloSitesElt->get_children(eltName::molInstanceRef);
      std::transform(instanceRefNodes.begin(),
		     instanceRefNodes.end(),
		     std::inserter(alloSiteMap,
				   alloSiteMap.begin()),
		     makeAlloSiteMapEntry(parsedPlex));

      // Install the allosteric sites into the plexFamily under the trivial
      // query which is always satisfied.  Later, when allosteric-plex
      // and allosteric-omni constructs include a state query, that query
      // will be installed here.
      pFamily->addAlloQueryAndMap(pPlexQueries,
				  alloSiteMap);

      //      pFamily->setAlloSites(alloSiteMap);
    }
  };

  // This function is for general parsing of a plex class (plex with state
  // specifications).  This could always have been used in parsing omni
  // dumpables and plex dumpables. It could now be used in in allosteric-plex
  // and allosteric-omni parsing, as well as in for the new bndKinase reaction
  // generator.
  //
  // It produces the recognized complex, together with a query which, when
  // applied to a species of the recognized complex, tells whether the species
  // satisfies the state specifications.
  plexFamily*
  parsePlexClass::operator()(xmlpp::Node* pParentNode) const
  {
    xmlpp::Element* pParentElt
      = domUtils::mustBeElementPtr(pParentNode);

    // Parse the complex.
    xmlpp::Element* pPlexElt
      = domUtils::mustGetUniqueChild(pParentElt,
				     eltName::plex);
    plexFamily* pFamily = recognizePlexElt(pPlexElt,
					   rPlex,
					   rMolUnit,
					   rPlexUnit);

    // Process the query specifications: a list of mol instances
    // with a "state filter" for each.  Now we have made this element
    // optional, so that allosteric-plex and allosteric-omni can use
    // it without their having to be changed.
    xmlpp::Node::NodeList instanceStatesNodes
      = pParentElt->get_children(eltName::instanceStates);

    if(instanceStatesNodes.size() != 0)
      {
	xmlpp::Element* pInstanceStatesElt
	  = domUtils::mustBeElementPtr(instanceStatesNodes.front());
	  
	// Eventually we will have to parse state filters for different
	// kinds of mols, but for now, just mod-mols
	xmlpp::Node::NodeList modMolInstanceRefNodes
	  = pInstanceStatesElt->get_children
	  (eltName::modMolInstanceRef);

	std::for_each(modMolInstanceRefNodes.begin(),
		      modMolInstanceRefNodes.end(),
		      addInstanceClausesToQuery(pQuery,
						rPlex,
						rMolUnit,
						rPlexUnit));
      }
    return pFamily;
  }

  class addPlexDumpableToFamily : public
  std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plexUnit& rPlexUnit;
    
  public:
    addPlexDumpableToFamily(mzr::mzrUnit& refMzrUnit,
			    bnd::molUnit& refMolUnit,
			    plexUnit& refPlexUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}
    
    void
    operator()(xmlpp::Node* pPlexSpeciesStreamNode) const throw(std::exception)
    {
      // Make sure that the node is an element, possibly unnecessarily.
      xmlpp::Element* pPlexSpeciesStreamElt
	= domUtils::mustBeElementPtr(pPlexSpeciesStreamNode);

      // Get the name of the species stream.
      std::string streamName
	= domUtils::mustGetAttrString
	(pPlexSpeciesStreamElt,
	 eltName::plexSpeciesStream_nameAttr);

      // Process the complex structure.
      parserPlex parsedPlex;
      xmlpp::Element* pPlexElt
	= domUtils::mustGetUniqueChild(pPlexSpeciesStreamElt,
				       eltName::plex);
      plexFamily* pFamily = recognizePlexElt(pPlexElt,
					     parsedPlex,
					     rMolUnit,
					     rPlexUnit);

      // Process the query specifications: a list of mol instances
      // with a "state filter" for each.  Now we have made this element
      // optional, so that allosteric-plex and allosteric-omni can use
      // it without their having to be changed.
      xmlpp::Node::NodeList instanceStatesNodes
	= pPlexSpeciesStreamElt->get_children(eltName::instanceStates);

      // Create the query that decides whether a plexSpecies's population
      // is to be included in this dumpable, and register it for
      // memory management by the plexUnit.
      andPlexQueries* pAndPlexQueries = new andPlexQueries();
      rPlexUnit.addPlexQuery(pAndPlexQueries);
      
      if(instanceStatesNodes.size() != 0)
	{
	  xmlpp::Element* pInstanceStatesElt
	    = domUtils::mustBeElementPtr(instanceStatesNodes.front());

	  // Eventually we will have to parse state filters for different
	  // kinds of mols, but for now, just mod-mols
	  xmlpp::Node::NodeList modMolInstanceRefNodes
	    = pInstanceStatesElt->get_children
	    (eltName::modMolInstanceRef);

	  // Parse mol instance state specifications into clauses
	  // in the "filter query" for the dumpable.  Plex species
	  // that pass the query test are added to the dumpable's
	  // list of species to dump.
	  std::for_each(modMolInstanceRefNodes.begin(),
			modMolInstanceRefNodes.end(),
			addInstanceClausesToQuery(pAndPlexQueries,
						  parsedPlex,
						  rMolUnit,
						  rPlexUnit));
	}

      // Create the dumpable.
      mzr::paramDumpable<plexSpecies>* pDumpable
	= new mzr::paramDumpable<plexSpecies>(streamName,
					      pAndPlexQueries);

      // Add the dumpable to the plexFamily.
      pFamily->addDumpable(pDumpable);

      // Add the dumpable to moleculizer for deletion and access.
      rMzrUnit.addSpeciesDumpable(pDumpable);
    }
  };

  // Transforms an XML modification map into a C++ version of the same
  // thing: a mapping from modification site names to modifications.
  //
  // The main action here is to look up the modification by name.
  class transformModification : public std::unary_function
  <const xmlpp::Node*, std::pair<std::string, const bnd::modification*> >
  {
    bnd::molUnit& rMolUnit;
    
  public:
    transformModification(bnd::molUnit& refMolUnit) :
      rMolUnit(refMolUnit)
    {}
    
    std::pair<std::string, const bnd::modification*>
    operator()(const xmlpp::Node* pModSiteRefNode) const throw(mzr::mzrXcpt)
    {
      // Cast the node pointer to an element pointer; probably unnecessarily
      // dynamically.
      const xmlpp::Element* pModSiteRefElt
	= domUtils::mustBeElementPtr(pModSiteRefNode);

      // Get the modification site name.
      std::string modSiteName
	= domUtils::mustGetAttrString(pModSiteRefElt,
				      bnd::eltName::modSiteRef_nameAttr);

      // Get the modification element, and from it, the name of the
      // modification.
      xmlpp::Element* pModRefElt
	= domUtils::mustGetUniqueChild(pModSiteRefElt,
				       bnd::eltName::modRef);
      std::string modName
	= domUtils::mustGetAttrString(pModRefElt,
				      bnd::eltName::modRef_nameAttr);

      // Look up the modification in the BAD BAD BAD static catalog
      // of modifications.
      return std::make_pair(modSiteName,
			    rMolUnit.getMod(modName));
    }
  };

  void
  replaceModMolDefaultState::operator()
    (const xmlpp::Node* pModMolInstanceRefNode) const throw(mzr::mzrXcpt)
  {
    // Make sure the mod-mol-instance-ref is an element; probably
    // unnecessarily dynamic cast.
    const xmlpp::Element* pModMolInstanceRefElt
      = dynamic_cast<const xmlpp::Element*>(pModMolInstanceRefNode);
    if(0 == pModMolInstanceRefElt)
      throw domUtils::badElementCastXcpt(pModMolInstanceRefNode);

    // Get the instance name of the mod-mol instance.
    std::string instanceName
      = domUtils::mustGetAttrString
      (pModMolInstanceRefElt,
       eltName::modMolInstanceRef_nameAttr);

    // Get the index of the mod-mol instance.
    int instanceNdx = rPlex.getMolNdxByName(instanceName);

    // Get the mod-mol itself and verify that it is a mod-mol.
    const plex& rParadigm = pFamily->getParadigm();
    bnd::modMol* pModMol
      = dynamic_cast<bnd::modMol*>(rParadigm.mols[instanceNdx]);
    if(0 == pModMol)
      throw instanceNotModMolXcpt(pModMolInstanceRefElt,
				  instanceName);

    // Target for parsing mod-map elements.  The mod mol
    // can intern this to to generate the corresponding modMolState*.
    std::map<std::string, const bnd::modification*> modMap;

    // Parse the modifcation map.
    xmlpp::Element* pModMapElt
      = domUtils::mustGetUniqueChild(pModMolInstanceRefElt,
				     bnd::eltName::modMap);
    xmlpp::Node::NodeList modSiteRefs
      = pModMapElt->get_children(bnd::eltName::modSiteRef);
    std::transform(modSiteRefs.begin(),
		   modSiteRefs.end(),
		   std::inserter(modMap,
				 modMap.begin()),
		   transformModification(rMolUnit));

    // Intern the modMolState and insert it into the vector of
    // molParams.  (This vector of molParams is the essential
    // stuff to construct the plexSpecies, already given its
    // structure.)
    rParams[instanceNdx] = pModMol->internModMap(modMap);
  }

  // Processes the combination of plex/instance-states elements to get
  // a plex species. (This combination can also be parsed to get a
  // plexFamily and an "andPlexQuery".)  This can be used to parse
  // explicit plex species (which also have a name) and plex species
  // in other contexts (for units that involve particular plex species
  // as substrates, etc.)
  namespace 
  {
    plexSpecies*
    plexSpeciesForClass(xmlpp::Element* pPlexElt,
			xmlpp::Node::NodeList& rModMolInstanceRefs,
			bnd::molUnit& rMolUnit,
			plexUnit& rPlexUnit)
      throw(std::exception)
    {
      // Get the plexFamily and arrange for remapping mol instance names
      // to indices in the paradigm plex.
      parserPlex thePlex;
      plexFamily* pFamily = recognizePlexElt(pPlexElt,
					     thePlex,
					     rMolUnit,
					     rPlexUnit);

      // Get vector of default states for the instance mols in the paradigm
      // as the first cut of the mol params of the plex species.
      std::vector<bnd::molParam> molParams
	= pFamily->makeDefaultMolParams();

      // Parse the various different kinds of state specification.
      //
      // For now, there are only mod-mols.
//       xmlpp::Node::NodeList modMolInstanceRefs
// 	= pInstanceStatesElt->get_children
// 	(eltName::modMolInstanceRef);

      std::for_each(rModMolInstanceRefs.begin(),
		    rModMolInstanceRefs.end(),
		    replaceModMolDefaultState(molParams,
					      pFamily,
					      thePlex,
					      rMolUnit));

      // Construct the plexSpecies.
      return pFamily->getMember(molParams);
    }
  }

  // Constructs the explicit plex species.  This has to happen after
  // all the omniplexes have been recognized, all the allosteric modifications
  // associated with plexes and omniplexes installed, and all the plexes
  // have been connected to their features.
  class processPlexSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plexUnit& rPlexUnit;
    
  public:
    processPlexSpecies(mzr::mzrUnit& refMzrUnit,
		       bnd::molUnit& refMolUnit,
		       plexUnit& refPlexUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}
    
    void
    operator()(xmlpp::Node* pPlexSpeciesNode) const throw(mzr::mzrXcpt)
    {
      // Make sure the plex-species node is an element; probably unnecessarily
      // dynamic.
      xmlpp::Element* pPlexSpeciesElt
	= dynamic_cast<xmlpp::Element*>(pPlexSpeciesNode);
      if(0 == pPlexSpeciesElt)
	throw domUtils::badElementCastXcpt(pPlexSpeciesNode);

      // Get the name of the plex species, for installation into the
      // BAD BAD BAD static catalog of named species.
      std::string plexSpeciesName
	= domUtils::mustGetAttrString(pPlexSpeciesElt,
				      eltName::plexSpecies_nameAttr);


      // Get the plexFamily and arrange for remapping mol instance names
      // to indices in the paradigm plex.
      xmlpp::Element* pPlexElt
	= domUtils::mustGetUniqueChild(pPlexSpeciesElt,
				       eltName::plex);
      parserPlex thePlex;
      recognizePlexElt(pPlexElt,
		       thePlex,
		       rMolUnit,
		       rPlexUnit);

      // Parse the states of the instance mols.  The syntax of this is the same
      // as the syntax of the instance state filter in a plex species stream,
      // but here, a modification site that is not mentioned is assumed to have
      // its default modification.  There (addInstanceClausesToQuery), an
      // unmentioned modification site was taken as a "wildcard", i.e. didn't
      // contribute to the filter.

      // Get the instance states, if any.  This is rather hacked up, since the
      // "instance-states" element is really useless, except for helping the
      // user with a little structure, and is now optional, due to the changes
      // in allosteric-plex and allosteric-omni.  This node list should be
      // empty when the "instance-states" element is not supplied.
      //
      // It's a hack because only modMol states can be done this way,
      // but there the only mols that support state right now.
      xmlpp::Node::NodeList modMolInstanceRefs;
      {
	xmlpp::Node::NodeList instanceStatesNodes
	  = pPlexSpeciesElt->get_children(eltName::instanceStates);
	if(instanceStatesNodes.size())
	  {
	    xmlpp::Element* pInstanceStatesElt
	      = domUtils::mustBeElementPtr(instanceStatesNodes.front());

	    modMolInstanceRefs
	      = pInstanceStatesElt->get_children(eltName::modMolInstanceRef);
	  }
      }
	  
      // Parse the particular plex species.
      plexSpecies* pPlexSpecies = plexSpeciesForClass(pPlexElt,
						      modMolInstanceRefs,
						      rMolUnit,
						      rPlexUnit);

      // Insert the plexSpecies into the BAD BAD BAD static catalog of named
      // species.
      rMzrUnit.addUserSpecies(plexSpeciesName,
			      pPlexSpecies);
    }
  };

  // Constructs the tagged plex species.  This has to happen after
  // all the omniplexes have been recognized, all the allosteric modifications
  // associated with plexes and omniplexes installed, and all the plexes
  // have been connected to their features.
  //
  // This is the same as the treatment of explicit plex species, sans the
  // name.
  class processTaggedPlexSpecies :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plexUnit& rPlexUnit;
    std::vector<plexSpecies*>& rUpdated;
    
  public:
    processTaggedPlexSpecies(mzr::mzrUnit& refMzrUnit,
			     bnd::molUnit& refMolUnit,
			     plexUnit& refPlexUnit,
			     std::vector<plexSpecies*>& rUpdatedSpecies) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      rUpdated(rUpdatedSpecies)
    {}
    
    void
    operator()(xmlpp::Node* pTaggedPlexSpeciesNode) const
      throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pTaggedPlexSpeciesElt
	= domUtils::mustBeElementPtr(pTaggedPlexSpeciesNode);

      // Get the plexFamily and arrange for remapping mol instance names
      // to indices in the paradigm plex.
      xmlpp::Element* pPlexElt
	= domUtils::mustGetUniqueChild(pTaggedPlexSpeciesElt,
				       eltName::plex);
      parserPlex thePlex;
      recognizePlexElt(pPlexElt,
		       thePlex,
		       rMolUnit,
		       rPlexUnit);

      // Parse the states of the instance mols.  The syntax of this is the same
      // as the syntax of the instance state filter in a plex species stream,
      // but here, a modification site that is not mentioned is assumed to have
      // its default modification.  There (addInstanceClausesToQuery), an
      // unmentioned modification site was taken as a "wildcard", i.e. didn't
      // contribute to the filter.

      // Get the instance states, if any.  This is rather hacked up, since the
      // "instance-states" element is really useless, except for helping the
      // user with a little structure, and is now optional, due to the changes
      // in allosteric-plex and allosteric-omni.  This node list should be
      // empty when the "instance-states" element is not supplied.
      //
      // It's a hack because only modMol states can be done this way,
      // but there the only mols that support state right now.
      xmlpp::Node::NodeList modMolInstanceRefs;
      {
	xmlpp::Node::NodeList instanceStatesNodes
	  = pTaggedPlexSpeciesElt->get_children(eltName::instanceStates);

	if(instanceStatesNodes.size())
	  {
	    xmlpp::Element* pInstanceStatesElt
	      = domUtils::mustBeElementPtr(instanceStatesNodes.front());

	    modMolInstanceRefs
	      = pInstanceStatesElt->get_children(eltName::modMolInstanceRef);
	  }
      }

      // Parse the particular plex species.  The species should cause
      // notifications during "prepareToRun" generating all their associated
      // reactions.  (We don't have any use of the species itself at this
      // point, unlike in processPlexSpecies.)
      plexSpecies* pSpecies
	= plexSpeciesForClass(pPlexElt,
			      modMolInstanceRefs,
			      rMolUnit,
			      rPlexUnit);

      // Had the plexSpecies ever been updated before the state dump?
      // If so, put it on the list for "zero-update."
      xmlpp::Node::NodeList updatedNodes
	= pTaggedPlexSpeciesElt->get_children(eltName::updated);
      if(0 < updatedNodes.size()) rUpdated.push_back(pSpecies);
      
    }
  };

  class processTaggedPlexSpeciesNpop :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plexUnit& rPlexUnit;
    std::set<mzr::reaction*>& rAffected;
    
  public:
    processTaggedPlexSpeciesNpop
    (mzr::mzrUnit& refMzrUnit,
     bnd::molUnit& refMolUnit,
     plexUnit& refPlexUnit,
     std::set<mzr::reaction*>& rAffectedReactions) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      rAffected(rAffectedReactions)
    {}
    
    void
    operator()(xmlpp::Node* pTaggedPlexSpeciesNode) const
      throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pTaggedPlexSpeciesElt
	= domUtils::mustBeElementPtr(pTaggedPlexSpeciesNode);

      // Get the plexFamily and arrange for remapping mol instance names
      // to indices in the paradigm plex.
      xmlpp::Element* pPlexElt
	= domUtils::mustGetUniqueChild(pTaggedPlexSpeciesElt,
				       eltName::plex);
      parserPlex thePlex;
      recognizePlexElt(pPlexElt,
		       thePlex,
		       rMolUnit,
		       rPlexUnit);

      // Parse the states of the instance mols.  The syntax of this is the same
      // as the syntax of the instance state filter in a plex species stream,
      // but here, a modification site that is not mentioned is assumed to have
      // its default modification.  There (addInstanceClausesToQuery), an
      // unmentioned modification site was taken as a "wildcard", i.e. didn't
      // contribute to the filter.

      // Get the instance states, if any.  This is rather hacked up, since the
      // "instance-states" element is really useless, except for helping the
      // user with a little structure, and is now optional, due to the changes
      // in allosteric-plex and allosteric-omni.  This node list should be
      // empty when the "instance-states" element is not supplied.
      //
      // It's a hack because only modMol states can be done this way,
      // but there the only mols that support state right now.
      xmlpp::Node::NodeList modMolInstanceRefs;
      {
	xmlpp::Node::NodeList instanceStatesNodes
	  = pTaggedPlexSpeciesElt->get_children(eltName::instanceStates);

	if(instanceStatesNodes.size())
	  {
	    xmlpp::Element* pInstanceStatesElt
	      = domUtils::mustBeElementPtr(instanceStatesNodes.front());

	    modMolInstanceRefs
	      = pInstanceStatesElt->get_children(eltName::modMolInstanceRef);
	  }
      }

      // Parse the particular plex species.  The species should cause
      // notifications during "prepareToRun" generating all their associated
      // reactions.  (We don't have any use of the species itself at this
      // point, unlike in processPlexSpecies.)
      plexSpecies* pSpecies
	= plexSpeciesForClass(pPlexElt,
			      modMolInstanceRefs,
			      rMolUnit,
			      rPlexUnit);

      // Get the population of the species at dump time.
      xmlpp::Element* pPopulationElt
	= domUtils::mustGetUniqueChild(pTaggedPlexSpeciesElt,
				       eltName::population);

      int dumpedPop
	= domUtils::mustGetAttrInt(pPopulationElt,
				   eltName::population_countAttr);

      // Update the species with its dumped population.  This obviates
      // the "zeroUpdate" pass that has to be made in prepareToDump,
      // as well as the usual "prepareToRun" which updates the
      // explicit plex species with the populations provided in
      // moleculizer-input.
      xmlpp::Node::NodeList updatedNodes
	= pTaggedPlexSpeciesElt->get_children(eltName::updated);
      if(0 < updatedNodes.size())
	{
	  pSpecies->update(rAffected,
			   dumpedPop);
	}
    }
  };

  // Does minimal recognition of an omniplex, without the usual initialization
  // of the resulting plexFamily.  This is for the preliminary scan to find
  // all the omniplexes.
  class unifyOmniplexNode : public
  std::unary_function<xmlpp::Node*, plexFamily*>
  {
    bnd::molUnit& rMolUnit;
    plexUnit& rPlexUnit;
    
  public:
    unifyOmniplexNode(bnd::molUnit& refMolUnit,
		      plexUnit& refPlexUnit) :
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}
    
    plexFamily*
    operator()(xmlpp::Node* pPlexNode) const throw(std::exception)
    {
      plexFamily* pFamily
	= unifyPlexElt(domUtils::mustBeElementPtr(pPlexNode),
		       rMolUnit,
		       rPlexUnit);

      return pFamily;
    }
  };

  // It seems to me that I should be able to do this with std::ptr_fun and
  // std::bind2nd, but there seems to be a problem with using a binder on a
  // function taking a reference argument.  (The binder tries to have a member
  // variable that is a reference to the argument type, which in this case is
  // a reference, and you can't make a reference to a reference.  This just
  // seems broken about STL.
  class mustGetUniqueChildWithName : public
  std::unary_function<xmlpp::Node*, xmlpp::Element*>
  {
    const std::string& rName;
  public:
    mustGetUniqueChildWithName(const std::string& rElementName) :
      rName(rElementName)
    {}
  
    xmlpp::Element*
    operator()(xmlpp::Node* pNode) const throw(std::exception)
    {
      return domUtils::mustGetUniqueChild(pNode,
					  rName);
    }
  };

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
    operator()(const std::string& rXpath) const throw(std::exception)
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

  void
  plexUnit::parseDomInput(xmlpp::Element* pRootElt,
			  xmlpp::Element* pModelElt,
			  xmlpp::Element* pStreamsElt,
			  xmlpp::Element* pEventsElt) throw(std::exception)
  {
    // First we pick out a number of elements and lists of elements
    // for tranversal.

    // Species streams.
    xmlpp::Element* pSpeciesStreamsElt
      = domUtils::mustGetUniqueChild(pStreamsElt,
				     mzr::eltName::speciesStreams);

    xmlpp::Node::NodeList omniSpeciesStreamNodes
      = pSpeciesStreamsElt
      ->get_children(eltName::omniSpeciesStream);

    xmlpp::Node::NodeList plexSpeciesStreamNodes
      = pSpeciesStreamsElt
      ->get_children(eltName::plexSpeciesStream);

    // Allostery specifications.
    // 
    // Is this element dispensable?
    xmlpp::Element* pAlloOmnisElt
      = domUtils::mustGetUniqueChild(pModelElt,
				     eltName::allostericOmnis);
    xmlpp::Node::NodeList alloOmniNodes
      = pAlloOmnisElt->get_children(eltName::allostericOmni);

    // Is this element dispensable?
    xmlpp::Element* pAlloPlexesElt
      = domUtils::mustGetUniqueChild(pModelElt,
				     eltName::allostericPlexes);
    xmlpp::Node::NodeList alloPlexNodes
      = pAlloPlexesElt->get_children(eltName::allostericPlex);

    // The explicit plexSpecies.
    xmlpp::Element* pExplicitSpeciesElt
      = domUtils::mustGetUniqueChild(pModelElt,
				     mzr::eltName::explicitSpecies);
    xmlpp::Node::NodeList plexSpeciesNodes
      = pExplicitSpeciesElt->get_children(eltName::plexSpecies);

    xmlpp::Node::NodeList omniPlexNodes;
    std::for_each(omniXpaths.begin(),
		  omniXpaths.end(),
		  addPathNodesToList(omniPlexNodes,
				     pRootElt));

    // "Unify" omniplex families (recognize, but without the usual
    // initializations) and put them on the plexUnit's list of omniplexes.
    // After this is done, plexes and omniplexes can be recognized in the usual
    // way.
    std::transform(omniPlexNodes.begin(),
		   omniPlexNodes.end(),
		   std::inserter(omniPlexFamilies,
				 omniPlexFamilies.begin()),
		   unifyOmniplexNode(rMolUnit,
				     *this));

    // Since the omniplex families have been "unified" in, they won't
    // undergo normal initialization when recognized.  Hence, we
    // run through them all and connect them to their features.
    std::for_each(omniPlexFamilies.begin(),
		  omniPlexFamilies.end(),
		  std::mem_fun(&plexFamily::connectToFeatures));

    // Do allosteric modifications of complexes, both omni and "regular".
    // These can be done with the same function, since allosteric-omnis
    // and allosteric-plexes have the same content, defined in the schema
    // by "allo-plex-conetent."
    // 
    // This must be completed before any species of complexes are generated.
    std::for_each(alloOmniNodes.begin(),
		  alloOmniNodes.end(),
		  parseAllostericOmni(rMolUnit,
				      *this));
				      
//     std::for_each(alloOmniNodes.begin(),
// 		  alloOmniNodes.end(),
// 		  doAlloMods(rMolUnit,
// 			     *this));
    std::for_each(alloPlexNodes.begin(),
		  alloPlexNodes.end(),
		  doAlloMods(rMolUnit,
			     *this));

    // Attach dumpables to families of complexes.  This must be done before any
    // species of complexes are generated.

    // Handling of query-based dumpables is slightly different for
    // omni-species-streams and plex-species streams.  (A plexFamily may have
    // both omni-species-streams and plex-species-streams attached to it, for
    // example.)

    // Parse query-based dumpables for omniplexes.
    std::for_each(omniSpeciesStreamNodes.begin(),
		  omniSpeciesStreamNodes.end(),
		  parseOmniSpeciesStream(rMzrUnit,
					 rMolUnit,
					 *this));

    // Parse query-based dumpables for plexes.
    std::for_each(plexSpeciesStreamNodes.begin(),
		  plexSpeciesStreamNodes.end(),
		  addPlexDumpableToFamily(rMzrUnit,
					  rMolUnit,
					  *this));

    // Note that all reaction generators should be in place before we generate
    // any species at all.  Features need to be notified of all species, since
    // features maintain lists of all the species that display the feature,
    // in order to generate reactions for new species.  This would imply
    // that any listing of reactions that we do in the state document will
    // be for other programs; Moleculizer won't be able to integrate them
    // in general.
    std::for_each(plexSpeciesNodes.begin(),
		  plexSpeciesNodes.end(),
		  processPlexSpecies(rMzrUnit,
				     rMolUnit,
				     *this));
  }

  class createInitialPop : public
  std::unary_function<xmlpp::Node*, void>
  {
    mzr::moleculizer& rMolzer;
    mzr::mzrUnit& rMzrUnit;
  public:
    createInitialPop(mzr::moleculizer& rMoleculizer,
		     mzr::mzrUnit& rMzr) :
      rMolzer(rMoleculizer),
      rMzrUnit(rMzr)
    {}
    
    void
    operator()(xmlpp::Node* pPlexSpeciesNode) const
      throw(std::exception)
    {
      xmlpp::Element* pPlexSpeciesElt
	= domUtils::mustBeElementPtr(pPlexSpeciesNode);

      std::string speciesName
	= domUtils::mustGetAttrString(pPlexSpeciesElt,
				      eltName::plexSpecies_nameAttr);

      xmlpp::Element* pPopulationElt
	= domUtils::mustGetUniqueChild(pPlexSpeciesElt,
				       mzr::eltName::population);
      int pop
	= domUtils::mustGetAttrInt(pPopulationElt,
				   mzr::eltName::population_countAttr);
      mzr::species* pSpecies
	= rMzrUnit.findSpecies(speciesName);

      if(0 == pSpecies)
	throw mzr::unknownSpeciesXcpt(pPlexSpeciesElt,
				      speciesName);

      // Rather than scheduling an event, we effectively do a create event.
      mzr::createEvent creator(pSpecies,
			       pop,
			       rMzrUnit);
      creator.doEvent(rMolzer,
		      rMzrUnit);
    }
  };

  void
  plexUnit::prepareToRun(xmlpp::Element* pRootElt,
			 xmlpp::Element* pModelElt,
			 xmlpp::Element* pStreamsElt,
			 xmlpp::Element* pEventsElt)
    throw(std::exception)
  {
    // No need for ensuring that all the plexSpecies have been
    // had their notification routines run.  This happens on
    // first update now.

    // Create the initial population of all explicit plex species.
    xmlpp::Element* pExplicitSpeciesElt
      = domUtils::mustGetUniqueChild(pModelElt,
				     mzr::eltName::explicitSpecies);
    xmlpp::Node::NodeList plexSpeciesNodes
      = pExplicitSpeciesElt->get_children(eltName::plexSpecies);

    // Since the species are explicit, we don't have to parse them
    // again, just look them up, to create their initial population.
    std::for_each(plexSpeciesNodes.begin(),
		  plexSpeciesNodes.end(),
		  createInitialPop(rMolzer,
				   rMzrUnit));
    
  }

//   class zeroUpdate : public
//   std::unary_function<std::multimap<int, plexFamily*>::value_type, void>
//   {
//     std::set<mzr::reaction*>& rAffectedReactions;
//   public:
//     zeroUpdate(std::set<mzr::reaction*>& rAffected) :
//       rAffectedReactions(rAffected)
//     {}

//     void
//     operator()(const argument_type& rIntFamilyPair) const
//     {
//       rIntFamilyPair.second->zeroUpdate(rAffectedReactions);
//     }
//   };

  class accumulateSpecies : public
  std::unary_function<std::multimap<int, plexFamily*>::value_type, void>
  {
    std::vector<plexSpecies*>& rAllSpecies;
  public:
    accumulateSpecies(std::vector<plexSpecies*>& rAllSpeciesVector) :
      rAllSpecies(rAllSpeciesVector)
    {
    }

    void
    operator()(const argument_type& rEntry) const
    {
      plexFamily* pPlexFamily = rEntry.second;
      pPlexFamily->accumulateSpecies(rAllSpecies);
    }
  };

  class zeroUpdateSpecies : public
  std::unary_function<plexSpecies*, void>
  {
    std::set<mzr::reaction*>& rAffected;
    
  public:
    zeroUpdateSpecies(std::set<mzr::reaction*>& rAffectedReactions) :
      rAffected(rAffectedReactions)
    {
    }

    void
    operator()(plexSpecies* pPlexSpecies) const
    {
      pPlexSpecies->update(rAffected,
			   0);
    }
  };
  
  void
  plexUnit::prepareToDump(xmlpp::Element* pRootElt,
			  xmlpp::Element* pModelElt,
			  xmlpp::Element* pStreamsElt,
			  xmlpp::Element* pEventsElt,
			  xmlpp::Element* pTaggedSpeciesElement)
    throw(std::exception)
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
    std::vector<plexSpecies*> updatedSpecies;
    std::for_each(taggedPlexSpeciesNodes.begin(),
		  taggedPlexSpeciesNodes.end(),
		  processTaggedPlexSpecies(rMzrUnit,
					   rMolUnit,
					   *this,
					   updatedSpecies));

    std::cerr << "Parsed "
	      << taggedPlexSpeciesNodes.size()
	      << " taggedPlexSpeciesNodes."
	      << std::endl
	      << "Of these, "
	      << updatedSpecies.size()
	      << " species had been updated."
	      << std::endl;

    // Run through all the plexSpecies, and update each by 0.
    // This should cause all the reactions for the species to come
    // into existence.  This may be a few more reactions than
    // were dumped, since some of the dumped species (with population 0)
    // may not ever have been updated.


    // Update each plexSpecies by zero, accumulating the
    // set of affected reactions.
    std::set<mzr::reaction*> affectedReactions;
    std::for_each(updatedSpecies.begin(),
		  updatedSpecies.end(),
		  zeroUpdateSpecies(affectedReactions));

    // Checking on just how bad the inflation of species and reactions is.
    std::vector<plexSpecies*> afterPlexSpecies;
    std::for_each(recognize.plexHasher.begin(),
		  recognize.plexHasher.end(),
		  accumulateSpecies(afterPlexSpecies));

    std::cerr << "After zeroUpdate, there are "
	      << afterPlexSpecies.size()
	      << " plex species."
	      << std::endl;

    // For consistency's sake, rescheduling the affected reactions.
    // Since the update was by 0, this isn't really necessary.
    std::for_each(affectedReactions.begin(),
		  affectedReactions.end(),
		  mzr::scheduleReaction(rMolzer.eventQ,
					rMzrUnit));

    // This is unit's default implementation of prepareToDump.
    prepareToRun(pRootElt,
		 pModelElt,
		 pStreamsElt,
		 pEventsElt);
  }

  void
  plexUnit::prepareToContinue(xmlpp::Element* pRootElt,
			      xmlpp::Element* pModelElt,
			      xmlpp::Element* pStreamsElt,
			      xmlpp::Element* pEventsElt,
			      std::map<std::string, std::string>& rTagToName,
			      xmlpp::Element* pTaggedSpeciesElement)
    throw(std::exception)
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
    std::set<mzr::reaction*> affectedReactions;
    std::for_each(taggedPlexSpeciesNodes.begin(),
		  taggedPlexSpeciesNodes.end(),
		  processTaggedPlexSpeciesNpop(rMzrUnit,
					       rMolUnit,
					       *this,
					       affectedReactions));

    // Reschedule the affected reactions.  The current simulation time
    // should have been set to the dump time by the mzrUnit.
    std::for_each(affectedReactions.begin(),
		  affectedReactions.end(),
		  mzr::scheduleReaction(rMolzer.eventQ,
					rMzrUnit));

    // In this version, do NOT run prepareToRun here, as that
    // sets the populations of the explicit species and reschedules
    // reactions.  May want rearchitecture here.
  }
}
