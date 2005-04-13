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

#include "mol/modMol.hh"
#include "mol/modQuery.hh"
#include "plex/plexFamily.hh"
#include "plex/parserPlex.hh"
#include "plex/parsePlex.hh"

namespace plx
{
  void
  parseMolInstance::
  operator()(xmlpp::Node* pMolInstNode) const
    throw(mzr::mzrXcpt)
  {
    // Make sure the node is an element, possibly unnecessarily.
    xmlpp::Element* pMolInstElt
      = dynamic_cast<xmlpp::Element*>(pMolInstNode);
    if(0 == pMolInstElt)
      throw domUtils::badElementCastXcpt(pMolInstNode);

    // Get the mol instance name.
    std::string molInstName
      = domUtils::mustGetAttrString(pMolInstElt,
				    eltName::molInstance_nameAttr);

    // Get the element telling which mol the instance is.
    xmlpp::Element* pMolRefElt
      = domUtils::mustGetUniqueChild(pMolInstElt,
				     eltName::molRef);

    // Get the mol name.
    std::string molName
      = domUtils::mustGetAttrString(pMolRefElt,
				    eltName::molRef_nameAttr);

    // Look up the mol in the molUnit.
    bnd::mol* pMol = rMolUnit.findMol(molName);
    if(0 == pMol) throw bnd::unknownMolXcpt(pMolRefElt,
					    molName);

    // Add the instance to the parserPlex (whose only function
    // beyond "plex" is to record the map from instance name
    // to instance index.
    rParsedPlex.addMolByName(molInstName,
			     pMol);
  }

  plexSiteSpec
  parseBindingPartner::
  operator()(xmlpp::Node* pMolInstRefNode) const
    throw(mzr::mzrXcpt)
  {
    // Make sure the node is an element, possibly unnecessarily.
    xmlpp::Element* pMolInstRefElt
      = dynamic_cast<xmlpp::Element*>(pMolInstRefNode);
    if(0 == pMolInstRefElt)
      throw domUtils::badElementCastXcpt(pMolInstRefNode);

    // Get the mol instance name.
    std::string molInstName
      = domUtils::mustGetAttrString
      (pMolInstRefElt,
       eltName::molInstanceRef_nameAttr);

    // Get the mol instance index; this is half of the returned
    // plexSiteSpec.
    int molNdx = rParsedPlex.getMolNdxByName(molInstName);
    if(0 > molNdx) throw unknownMolInstanceXcpt(pMolInstRefElt,
						molInstName);

    // In order to look up the binding site's index, get the
    // named mol.
    bnd::mol* pMol = rParsedPlex.mols[molNdx];

    // Get the element naming the binding site on the above
    // mol instance.
    xmlpp::Element* pBindingSiteRefElt
      = domUtils::mustGetUniqueChild(pMolInstRefElt,
				     bnd::eltName::bindingSiteRef);
    std::string bindingSiteName
      = domUtils::mustGetAttrString
      (pBindingSiteRefElt,
       bnd::eltName::bindingSiteRef_nameAttr);

    // Get the binding site index from the mol.
    int bindingSiteNdx = -1;
    if(! pMol->findSite(bindingSiteName,
			bindingSiteNdx))
      throw bnd::unknownSiteXcpt(pBindingSiteRefElt,
				 bindingSiteName);

    // Construct the return value.
    plexSiteSpec returnValue(molNdx,
			     bindingSiteNdx);

    // Test that the site is not currently bound, and enter the site
    // (somewhat prematurely) into the set of bound sites.
    std::pair<std::set<plexSiteSpec>::iterator, bool> insertResult
      = rBoundSites.insert(returnValue);
    if(! insertResult.second)
      throw multiplyBoundSiteXcpt(pBindingSiteRefElt,
				  bindingSiteName,
				  molInstName);
    return returnValue;
  }

  void
  parseBinding::
  operator()(xmlpp::Node* pBindingNode) const
    throw(mzr::mzrXcpt)
  {
    // Make sure the node is an element, possibly unnecessarily.
    xmlpp::Element* pBindingElt
      = dynamic_cast<xmlpp::Element*>(pBindingNode);
    if(0 == pBindingElt) throw domUtils::badElementCastXcpt(pBindingNode);

    // Parse the two binding partners, both elements named
    // "mol-instance-ref".
    xmlpp::Node::NodeList partners
      = pBindingElt->get_children(eltName::molInstanceRef);
    if(2 != partners.size())
      throw domUtils::badChildCountXcpt(pBindingElt,
					eltName::molInstanceRef,
					2,
					partners.size());

    // Transform the list of two mol-instance-ref elements into
    // a vector of two plexSiteSpecs.  A plexSiteSpec gives a mol
    // instance index, and on the corresponding mol, a binding
    // site index.
    //
    // This transformation also tests the binding site to ensure
    // that it is not already bound, and inserts the binding
    // site (somewhat prematurely) into the set of bound sites.
    std::vector<plexSiteSpec> siteSpecs;
    std::transform(partners.begin(),
		   partners.end(),
		   std::back_inserter(siteSpecs),
		   parseBindingPartner(rParsedPlex,
				       rBoundSites));

    // Add the new binding to the plex.
    rParsedPlex.bindings.push_back(plexBinding(siteSpecs[0],
					 siteSpecs[1]));
  }

  void
  parsePlex::
  operator()(xmlpp::Element* pPlexElt) const
    throw(mzr::mzrXcpt)
  {
    // Process the mol instances.
    xmlpp::Node::NodeList molInstances
      = pPlexElt->get_children(eltName::molInstance);
    std::for_each(molInstances.begin(),
		  molInstances.end(),
		  parseMolInstance(rMolUnit,
				   rParsedPlex));

    // Process the bindings.
    //
    // This is to make sure that no site is in two bindings.
    std::set<plexSiteSpec> boundSites;
    xmlpp::Node::NodeList bindings
      = pPlexElt->get_children(eltName::binding);
    std::for_each(bindings.begin(),
		  bindings.end(),
		  parseBinding(rParsedPlex,
			       boundSites));
  }

  namespace
  {
    // Adjusts the instance indices in an instance-name to instance-index map
    // using recognition isomorphism.
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
  }

  plexFamily*
  unifyPlexNode(xmlpp::Node* pPlexNode,
		bnd::molUnit& rMolUnit,
		plexUnit& rPlexUnit,
		parserPlex& rParsedPlex)
    throw(mzr::mzrXcpt)
  {
    xmlpp::Element* pPlexElt =
      domUtils::mustBeElementPtr(pPlexNode);

    // Parse the plex.
    parsePlex plexParser(rMolUnit,
			 rParsedPlex);
    plexParser(pPlexElt);

    // Recognize the plex, but without doing any of the usual
    // initializations of the resulting plexFamily.
    plexFamily* pFamily = 0;
    plexIsoPair recognitionIsos;
    rPlexUnit.recognize.unify(rParsedPlex,
			      pFamily,
			      & recognitionIsos);

    // Reindex the instance-name to instance-index map.
    std::map<std::string, int> nameToMolNdx;
    std::transform(rParsedPlex.nameToMolNdx.begin(),
		   rParsedPlex.nameToMolNdx.end(),
		   std::inserter(nameToMolNdx,
				 nameToMolNdx.begin()),
		   remapMolNdx(recognitionIsos));

    // Rearrange the instances of the parsed plex to coincide
    // with the paradigm of the plexFamily.
    rParsedPlex = parserPlex(nameToMolNdx,
			     pFamily->getParadigm());

    return pFamily;
  }

  plexFamily*
  recognizePlexElt(xmlpp::Element* pPlexElt,
		   parserPlex& rParsedPlex,
		   bnd::molUnit& rMolUnit,
		   plexUnit& rPlexUnit)
  {
    // Parse the plex.
    parserPlex thePlex;
    parsePlex plexParser(rMolUnit,
			 thePlex);
    plexParser(pPlexElt);

    // Recognize the plex, but without doing any of the usual
    // initializations of the resulting plexFamily.
    plexIsoPair recognitionIsos;
    plexFamily* pFamily = rPlexUnit.recognize(thePlex,
					      recognitionIsos);

    // Reindex the instance-name to instance-index map.
    std::map<std::string, int> nameToMolNdx;
    std::transform(thePlex.nameToMolNdx.begin(),
		   thePlex.nameToMolNdx.end(),
		   std::inserter(nameToMolNdx,
				 nameToMolNdx.begin()),
		   remapMolNdx(recognitionIsos));

    rParsedPlex = parserPlex(nameToMolNdx,
			     pFamily->getParadigm());

    return pFamily;
  }

  namespace
  {
    // Parses a single mod-site-ref/modification pair, constructs
    // a query from it, and adds the query to a given, overall query
    // for a plex.
    class parseModSiteRefQuery :
      public std::unary_function<xmlpp::Node*, void>
    {
      bnd::molUnit& rMolUnit;
      plexUnit& rPlexUnit;
      andPlexQueries* pQuery;
      bnd::modMol* pMol;
      plexMolSpec instanceNdx;
    public:
      parseModSiteRefQuery(bnd::molUnit& refMolUnit,
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
      operator()(xmlpp::Node* pModSiteRefNode) const throw(mzr::mzrXcpt)
      {
	xmlpp::Element* pModSiteRefElt
	  = domUtils::mustBeElementPtr(pModSiteRefNode);
    
	// Get the modification site name.
	std::string modSiteName
	  = domUtils::mustGetAttrString(pModSiteRefElt,
					bnd::eltName::modSiteRef_nameAttr);

	// Look up the modification site on the mol, getting the
	// index of the modification site.
	int modSiteNdx = pMol->mustGetModSiteNdx(pModSiteRefElt,
						 modSiteName);

	// Get the modification name.
	xmlpp::Element* pModRefElt
	  = domUtils::mustGetUniqueChild(pModSiteRefNode,
					 bnd::eltName::modRef);
	std::string modName
	  = domUtils::mustGetAttrString(pModRefElt,
					bnd::eltName::modRef_nameAttr);

	// Look up the modification.
	const bnd::modification* pMod
	  = rMolUnit.mustGetMod(pModRefElt,
				modName);
	
	// Construct the query and register for memory management
	// by the plexUnit.
	bnd::modMixinStateQuery* pModQuery
	  = new bnd::modMixinStateQuery(rPlexUnit,
					rMolUnit,
					instanceNdx,
					modSiteNdx,
					pMod);

	// Add the query to the conjunction of queries to be applied to
	// species of complexes appearing in this stream.
	pQuery->addQuery(pModQuery);
      }
    };

    // Parses a single mod-mol-instance-ref element, adding the
    // resulting query to a given, overall query for a plex.
    class parseModMolInstanceStateQuery :
      public std::unary_function<xmlpp::Node*, void>
    {
      andPlexQueries* pQuery;
      const parserPlex& rPlex;
      bnd::molUnit& rMolUnit;
      plexUnit& rPlexUnit;
    
    public:
      parseModMolInstanceStateQuery(andPlexQueries* pAndPlexQueries,
				       const parserPlex& rParsedPlex,
				       bnd::molUnit& refMolUnit,
				       plexUnit& refPlexUnit) :
	pQuery(pAndPlexQueries),
	rPlex(rParsedPlex),
	rMolUnit(refMolUnit),
	rPlexUnit(refPlexUnit)
      {}

      void
      operator()(xmlpp::Node* pModMolInstanceRefNode) const
	throw(mzr::mzrXcpt)
      {
	xmlpp::Element* pModMolInstanceRefElt
	  = domUtils::mustBeElementPtr(pModMolInstanceRefNode);

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

	// Make sure the mol is a mod-mol.
	bnd::modMol* pModMol
	  = bnd::mustBeModMolPtr(pModMolInstanceRefElt,
				 rPlex.mols[instanceNdx]);

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
		      parseModSiteRefQuery(rMolUnit,
					   rPlexUnit,
					   pQuery,
					   pModMol,
					   instanceNdx));
      }
    };

    // Parses state queries for mod-mol instances in a plex, adding them
    // to an existing, overall state query for the plex.
    //
    // Note that other (future) kinds of mols could have other kinds of state
    // queries, which would have the same parent node (pInstanceStatesNode).
    // This routine picks out just the modMolInstanceRef's, which go with
    // mod-mols.
    void
    parseModMolInstanceStateQueries(xmlpp::Node* pInstanceStatesNode,
				    andPlexQueries* pQuery,
				    const parserPlex& rPlex,
				    bnd::molUnit& rMolUnit,
				    plexUnit& rPlexUnit)
      throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pInstanceStatesElt
	= domUtils::mustBeElementPtr(pInstanceStatesNode);
	  
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
		    parseModMolInstanceStateQuery(pQuery,
						  rPlex,
						  rMolUnit,
						  rPlexUnit));

    }
  }

  // This routine is really just the place where new parsers for
  // different kinds of mol state specifications could be installed.
  void
  parseInstanceStateQueries(xmlpp::Node* pInstanceStatesNode,
			    andPlexQueries* pQuery,
			    const parserPlex& rParserPlex,
			    bnd::molUnit& rMolUnit,
			    plexUnit& rPlexUnit)
    throw(mzr::mzrXcpt)
  {
    // Parse mod-mol state queries, adding them to the overall plex query.
    parseModMolInstanceStateQueries(pInstanceStatesNode,
				    pQuery,
				    rParserPlex,
				    rMolUnit,
				    rPlexUnit);

    // This is where parsers for new kinds of mol state queries could
    // be installed.
  }

  plexFamily*
  parsePlexClass::
  operator()(xmlpp::Node* pParentNode) const
    throw(mzr::mzrXcpt)
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
	parseInstanceStateQueries(instanceStatesNodes.front(),
				  pQuery,
				  rPlex,
				  rMolUnit,
				  rPlexUnit);
      }

    return pFamily;
  }

  namespace
  {
    // Parses a single mod-site-ref/modification pair, converting it
    // into an entry in a map from modification site name to modification.
    class parseModSiteRef :
      public std::unary_function
    <xmlpp::Node*, std::pair<std::string, const bnd::modification*> >
    {
      bnd::molUnit& rMolUnit;
    
    public:
      parseModSiteRef(bnd::molUnit& refMolUnit) :
	rMolUnit(refMolUnit)
      {}

      std::pair<std::string, const bnd::modification*>
      operator()(xmlpp::Node* pModSiteRefNode) const
	throw(mzr::mzrXcpt)
      {
	xmlpp::Element* pModSiteRefElt
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

	// Look up the modification in bnd::molUnit.
	return std::make_pair(modSiteName,
			      rMolUnit.mustGetMod(modName));
      }
    };

    class parseModMolInstanceState :
      public std::unary_function<xmlpp::Node*, void>
    {
      std::vector<bnd::molParam>& rParams;
      const parserPlex& rPlex;
      const plexFamily* pFamily;
      bnd::molUnit& rMolUnit;

    public:
      parseModMolInstanceState(std::vector<bnd::molParam>& rMolParams,
			       const parserPlex& rParserPlex,
			       const plexFamily* pPlexFamily,
			       bnd::molUnit& refMolUnit) :
	rParams(rMolParams),
	rPlex(rParserPlex),
	pFamily(pPlexFamily),
	rMolUnit(refMolUnit)
      {}

      void
      operator()(xmlpp::Node* pModMolInstanceRefNode) const
	throw(mzr::mzrXcpt)
      {
	// Make sure the mod-mol-instance-ref is an element; probably
	// unnecessarily dynamic cast.
	xmlpp::Element* pModMolInstanceRefElt
	  = domUtils::mustBeElementPtr(pModMolInstanceRefNode);

	// Get the instance name of the mod-mol instance.
	std::string instanceName
	  = domUtils::mustGetAttrString(pModMolInstanceRefElt,
					eltName::modMolInstanceRef_nameAttr);

	// Get the index of the mod-mol instance.
	int instanceNdx
	  = rPlex.getMolNdxByName(instanceName);

	// Get the mod-mol itself and verify that it is a mod-mol.
	const plex& rParadigm = pFamily->getParadigm();
	bnd::modMol* pModMol
	  = bnd::mustBeModMolPtr(pModMolInstanceRefElt,
				 rParadigm.mols[instanceNdx]);

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
		       parseModSiteRef(rMolUnit));

	// Intern the modMolState and insert it into the vector of
	// molParams.  (This vector of molParams is the essential
	// stuff to construct the plexSpecies, already given its
	// structure.)
	//
	// Perhaps mod-mols should support a similar function for
	// queries.
	rParams[instanceNdx] = pModMol->internModMap(modMap);
      }
    };
    
    void
    parseModMolInstanceStates(xmlpp::Node* pInstanceStatesNode,
			      std::vector<bnd::molParam>& rMolParams,
			      const parserPlex& rParserPlex,
			      const plexFamily* pPlexFamily,
			      bnd::molUnit& rMolUnit)
    {
      xmlpp::Element* pInstanceStatesElt
	= domUtils::mustBeElementPtr(pInstanceStatesNode);
	  
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
		    parseModMolInstanceState(rMolParams,
					     rParserPlex,
					     pPlexFamily,
					     rMolUnit));
    }
  }

  void
  parseInstanceStates(xmlpp::Node* pInstanceStatesNode,
		      std::vector<bnd::molParam>& rMolParams,
		      const parserPlex& rParserPlex,
		      const plexFamily* pPlexFamily,
		      bnd::molUnit& rMolUnit)
    throw(mzr::mzrXcpt)
  {
    // Parse mod-mol states, replacing the states in rMolParams.
    parseModMolInstanceStates(pInstanceStatesNode,
			      rMolParams,
			      rParserPlex,
			      pPlexFamily,
			      rMolUnit);

    // This is where new parsers for new kinds of mol states for new kinds
    // of mols could be installed.
  }

  plexSpecies*
  parsePlexSpecies(xmlpp::Element* pParentElt,
		   mzr::mzrUnit& rMzrUnit,
		   bnd::molUnit& rMolUnit,
		   plexUnit& rPlexUnit)
    throw(mzr::mzrXcpt)
  {
    // Get the plexFamily and arrange for remapping mol instance names
    // to indices in the paradigm plex.
    xmlpp::Element* pPlexElt
      = domUtils::mustGetUniqueChild(pParentElt,
				     eltName::plex);
    parserPlex thePlex;
    plexFamily* pPlexFamily
      = recognizePlexElt(pPlexElt,
			 thePlex,
			 rMolUnit,
			 rPlexUnit);

    // If no instance states are specified, then we use the default molParams.
    std::vector<bnd::molParam> molParams
      = pPlexFamily->makeDefaultMolParams();

    // Parse optional instanceStates of the instance mols. These will replace
    // some of the molParams.
    xmlpp::Element* pInstanceStatesElt
      = domUtils::getOptionalChild(pParentElt,
				   eltName::instanceStates);
    if(pInstanceStatesElt)
      {
	parseInstanceStates(pInstanceStatesElt,
			    molParams,
			    thePlex,
			    pPlexFamily,
			    rMolUnit);
      }
	  
    // Construct the new plexSpecies from the molParams using the
    // parsed plex's family.
    return pPlexFamily->getMember(molParams);
  }

  plexSpecies*
  parseExplicitPlexSpecies::
  operator()(xmlpp::Node* pPlexSpeciesNode) const
    throw(mzr::mzrXcpt)
  {
    xmlpp::Element* pPlexSpeciesElt
      = domUtils::mustBeElementPtr(pPlexSpeciesNode);

    // Get the name of the plex species.
    std::string plexSpeciesName
      = domUtils::mustGetAttrString(pPlexSpeciesElt,
				    eltName::plexSpecies_nameAttr);

    // Parse the plexSpecies.
    plexSpecies* pPlexSpecies
      = parsePlexSpecies(pPlexSpeciesElt,
			 rMzrUnit,
			 rMolUnit,
			 rPlexUnit);

    // Insert the plexSpecies into mzrUnit's catalog of named species.
    rMzrUnit.addUserSpecies(plexSpeciesName,
			    pPlexSpecies);

    return pPlexSpecies;
  }

  plexSpecies*
  parseTaggedPlexSpecies::
  operator()(xmlpp::Node* pTaggedPlexSpeciesNode) const
    throw(mzr::mzrXcpt)
  {
    xmlpp::Element* pTaggedPlexSpeciesElt
      = domUtils::mustBeElementPtr(pTaggedPlexSpeciesNode);

    // Parse the plexSpecies.
    plexSpecies* pPlexSpecies
      = parsePlexSpecies(pTaggedPlexSpeciesElt,
			 rMzrUnit,
			 rMolUnit,
			 rPlexUnit);

    // Had the plexSpecies ever been updated before state was dumped?
    xmlpp::Node::NodeList updatedNodes
      = pTaggedPlexSpeciesElt->get_children(eltName::updated);
    if(0 < updatedNodes.size()) rUpdated.push_back(pPlexSpecies);

    return pPlexSpecies;
  }

  namespace
  {
    // Parses a mol-instance-ref node in an allosteric-sites element
    // to make an entry for a siteToShapeMap.
    class parseMolInstanceRef :
      public std::unary_function<xmlpp::Node*, siteToShapeMap::value_type>
    {
      const parserPlex& rPlex;

    public:
      parseMolInstanceRef(const parserPlex& rParserPlex) :
	rPlex(rParserPlex)
      {}

      result_type
      operator()(argument_type pMolInstanceRefNode) const
	throw(mzr::mzrXcpt)
      {
	// Check that the node is an element node, possibly unnecessarily.
	xmlpp::Element* pMolInstanceRefElt
	  = domUtils::mustBeElementPtr(pMolInstanceRefNode);
    
	// Get the mol instance name.
	std::string molInstanceName
	  = domUtils::mustGetAttrString
	  (pMolInstanceRefElt,
	   eltName::molInstanceRef_nameAttr);

	// Get the corresponding mol index, part of the plexSiteSpec, in order
	// to make an entry in the allostery map of the plex definition.
	int molNdx
	  = rPlex.mustGetMolNdxByName(pMolInstanceRefElt,
				      molInstanceName);

	// Get the corresponding mol, in order to look up binding site by
	// name.
	bnd::mol* pMol = rPlex.mols[molNdx];

	// Get the child element that gives the site name.
	xmlpp::Element* pBindingSiteRefElt
	  = domUtils::mustGetUniqueChild(pMolInstanceRefElt,
					 bnd::eltName::bindingSiteRef);

	// Get the name of the binding site.
	std::string bindingSiteName
	  = domUtils::mustGetAttrString
	  (pBindingSiteRefElt,
	   bnd::eltName::bindingSiteRef_nameAttr);

	// Get the index of the binding site, part of the plexSiteSpec.
	int bindingSiteNdx
	  = pMol->mustGetSiteNdxByName(pBindingSiteRefElt,
				       bindingSiteName);

	// Get the actual binding site, to look up the shape.
	bnd::bindingSite& rBindingSite
	  = pMol->getSiteByNdx(bindingSiteNdx);
    
	// Get the child element that gives the shape of the binding
	// site when in this subcomplex.
	xmlpp::Element* pSiteShapeRefElt
	  = domUtils::mustGetUniqueChild(pBindingSiteRefElt,
					 bnd::eltName::siteShapeRef);

	// Get the name of the site shape.
	std::string siteShapeName
	  = domUtils::mustGetAttrString(pSiteShapeRefElt,
					bnd::eltName::siteShapeRef_nameAttr);

	// Look up the site shape.
	const bnd::siteShape* pSiteShape
	  = rBindingSite.mustGetParam(pSiteShapeRefElt,
				      pMol,
				      siteShapeName);

	// Return entry for allostery map.
	return std::make_pair(plexSiteSpec(molNdx,
					   bindingSiteNdx),
			      pSiteShape);
      }
    };
  }

  void
  parseAllostericSites::
  operator()(xmlpp::Element* pAlloSitesElt) const
    throw(mzr::mzrXcpt)
  {
    xmlpp::Node::NodeList molInstanceRefNodes
      = pAlloSitesElt->get_children(eltName::molInstanceRef);

    std::transform(molInstanceRefNodes.begin(),
		   molInstanceRefNodes.end(),
		   std::inserter(rSiteToShapeMap,
				 rSiteToShapeMap.begin()),
		   parseMolInstanceRef(rParserPlex));
  }

  void
  parseAllostericPlex::
  operator()(xmlpp::Node* pAlloPlexNode) const
    throw(mzr::mzrXcpt)
  {
    // Parse the plex class.
    parserPlex parsedPlex;
    andPlexQueries* pAndPlexQueries = new andPlexQueries(rPlexUnit);
    parsePlexClass classParser(rMolUnit,
			       rPlexUnit,
			       parsedPlex,
			       pAndPlexQueries);
    plexFamily* pFamily
      = classParser(pAlloPlexNode);

    // Parse the allosteric sites.
    siteToShapeMap alloSiteMap;
    parseAllostericSites siteParser(parsedPlex,
				    alloSiteMap);

    // Add the query and alloSiteMap to the plexFamily's alloStateList, which
    // is used to test new plexSpecies as they appear, applying allosteric
    // modifications to those that answer the query.
    pFamily->addAlloQueryAndMap(pAndPlexQueries,
				alloSiteMap);
  }

  void
  parsePlexSpeciesStream::
  operator()(xmlpp::Node* pPlexSpeciesStreamNode) const
    throw(mzr::mzrXcpt)
  {
    // Make sure that the node is an element, possibly unnecessarily.
    xmlpp::Element* pPlexSpeciesStreamElt
      = domUtils::mustBeElementPtr(pPlexSpeciesStreamNode);

    // Get the name of the species stream.
    std::string streamName
      = domUtils::mustGetAttrString
      (pPlexSpeciesStreamElt,
       eltName::plexSpeciesStream_nameAttr);

    // Parse the plex class, to get a plex, plexFamily, and a state query.
    parserPlex parsedPlex;
    andPlexQueries* pAndPlexQueries = new andPlexQueries(rPlexUnit);
    // Is this ever used in an STL algorithm? Maybe just a function?  Why
    // do I bother?
    parsePlexClass classParser(rMolUnit,
			       rPlexUnit,
			       parsedPlex,
			       pAndPlexQueries);
    plexFamily* pFamily = classParser(pPlexSpeciesStreamElt);

    // Create the dumpable.
    mzr::paramDumpable<plexSpecies>* pDumpable
      = new mzr::paramDumpable<plexSpecies>(streamName,
					    pAndPlexQueries);

    // Add the dumpable to the plexFamily.
    pFamily->addDumpable(pDumpable);

    // Add the dumpable to moleculizer for deletion and access.
    rMzrUnit.addSpeciesDumpable(pDumpable);
  }
}
