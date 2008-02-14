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

#include "utl/badElementCastXcpt.hh"
#include "utl/badChildCountXcpt.hh"
#include "cpx/modMolStateQuery.hh"
#include "cml/cmlEltName.hh"
#include "cml/unkMolXcpt.hh"
#include "cml/unkSiteXcpt.hh"
#include "cml/badModMolXcpt.hh"
#include "clx/cptPlexFamily.hh"
#include "clx/parserPlex.hh"
#include "clx/parsePlex.hh"
#include "clx/multBoundSiteXcpt.hh"
#include "clx/parsedPlexNotConnXcpt.hh"
#include "clx/cptPlexQueries.hh"

namespace clx
{
  void
  parseMolInstance::
  operator()(xmlpp::Node* pMolInstNode) const
    throw(utl::xcpt)
  {
    // Make sure the node is an element, possibly unnecessarily.
    xmlpp::Element* pMolInstElt
      = dynamic_cast<xmlpp::Element*>(pMolInstNode);
    if(0 == pMolInstElt)
      throw utl::dom::badElementCastXcpt(pMolInstNode);

    // Get the mol instance name.
    std::string molInstName
      = utl::dom::mustGetAttrString(pMolInstElt,
				    eltName::molInstance_nameAttr);

    // Get the element telling which mol the instance is.
    xmlpp::Element* pMolRefElt
      = utl::dom::mustGetUniqueChild(pMolInstElt,
				     eltName::molRef);

    // Get the mol name.
    std::string molName
      = utl::dom::mustGetAttrString(pMolRefElt,
				    eltName::molRef_nameAttr);

    // Look up the mol in the cmlUnit.
    cml::cptMol* pMol = rCmlUnit.findMol(molName);
    if(0 == pMol) throw cml::unkMolXcpt(molName,
					pMolRefElt);

    // Add the instance to the parserPlex (whose only function
    // beyond "plex" is to record the map from instance name
    // to instance index.
    rParsedPlex.addMolByName(molInstName,
			     pMol);
  }

  cpx::siteSpec
  parseBindingPartner::
  operator()(xmlpp::Node* pMolInstRefNode) const
    throw(utl::xcpt)
  {
    // Make sure the node is an element, possibly unnecessarily.
    xmlpp::Element* pMolInstRefElt
      = dynamic_cast<xmlpp::Element*>(pMolInstRefNode);
    if(0 == pMolInstRefElt)
      throw utl::dom::badElementCastXcpt(pMolInstRefNode);

    // Get the mol instance name.
    std::string molInstName
      = utl::dom::mustGetAttrString
      (pMolInstRefElt,
       eltName::molInstanceRef_nameAttr);

    // Get the mol instance index; this is half of the returned
    // cpx::siteSpec.
    int molNdx = rParsedPlex.getMolNdxByName(molInstName);
    if(0 > molNdx) throw unkMolInstXcpt(molInstName,
					pMolInstRefElt);

    // In order to look up the binding site's index, get the
    // named mol.
    cml::cptMol* pMol = rParsedPlex.mols[molNdx];

    // Get the element naming the binding site on the above
    // mol instance.
    xmlpp::Element* pBindingSiteRefElt
      = utl::dom::mustGetUniqueChild(pMolInstRefElt,
				     cml::eltName::bindingSiteRef);
    std::string bindingSiteName
      = utl::dom::mustGetAttrString
      (pBindingSiteRefElt,
       cml::eltName::bindingSiteRef_nameAttr);

    // Get the binding site index from the mol.
    int bindingSiteNdx = -1;
    if(! pMol->findSite(bindingSiteName,
			bindingSiteNdx))
      throw cml::unkSiteXcpt(bindingSiteName,
			     pBindingSiteRefElt);

    // Construct the return value.
    cpx::siteSpec returnValue(molNdx,
			     bindingSiteNdx);

    // Test that the site is not currently bound, and enter the site
    // (somewhat prematurely) into the set of bound sites.
    std::pair<std::set<cpx::siteSpec>::iterator, bool> insertResult
      = rBoundSites.insert(returnValue);
    if(! insertResult.second)
      throw multBoundSiteXcpt(bindingSiteName,
			      molInstName,
			      pBindingSiteRefElt);
    return returnValue;
  }

  void
  parseBinding::
  operator()(xmlpp::Node* pBindingNode) const
    throw(utl::xcpt)
  {
    // Make sure the node is an element, possibly unnecessarily.
    xmlpp::Element* pBindingElt
      = dynamic_cast<xmlpp::Element*>(pBindingNode);
    if(0 == pBindingElt) throw utl::dom::badElementCastXcpt(pBindingNode);

    // Parse the two binding partners, both elements named
    // "mol-instance-ref".
    xmlpp::Node::NodeList partners
      = pBindingElt->get_children(eltName::molInstanceRef);
    if(2 != partners.size())
      throw utl::dom::badChildCountXcpt::general(pBindingElt,
						 eltName::molInstanceRef,
						 2,
						 partners.size());

    // Transform the list of two mol-instance-ref elements into
    // a vector of two cpx::siteSpecs.  A cpx::siteSpec gives a mol
    // instance index, and on the corresponding mol, a binding
    // site index.
    //
    // This transformation also tests the binding site to ensure
    // that it is not already bound, and inserts the binding
    // site (somewhat prematurely) into the set of bound sites.
    std::vector<cpx::siteSpec> siteSpecs;
    std::transform(partners.begin(),
		   partners.end(),
		   std::back_inserter(siteSpecs),
		   parseBindingPartner(rParsedPlex,
				       rBoundSites));

    // Add the new binding to the plex.
    rParsedPlex.bindings.push_back(cpx::binding(siteSpecs[0],
						siteSpecs[1]));
  }

  void
  parsePlex::
  operator()(xmlpp::Element* pPlexElt) const
    throw(utl::xcpt)
  {
    // Process the mol instances.
    xmlpp::Node::NodeList molInstances
      = pPlexElt->get_children(eltName::molInstance);
    std::for_each(molInstances.begin(),
		  molInstances.end(),
		  parseMolInstance(rCmlUnit,
				   rParsedPlex));

    // Process the bindings.
    //
    // This is to make sure that no site is in two bindings.
    std::set<cpx::siteSpec> boundSites;
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
      remapMolNdx(const cpx::plexIso& rIsoPair) :
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

  cptPlexFamily*
  unifyPlexNode(xmlpp::Node* pPlexNode,
		cml::cmlUnit& rCmlUnit,
		clxUnit& rClxUnit,
		parserPlex& rParsedPlex)
    throw(utl::xcpt)
  {
    xmlpp::Element* pPlexElt =
      utl::dom::mustBeElementPtr(pPlexNode);

    // Parse the plex.
    parsePlex plexParser(rCmlUnit,
			 rParsedPlex);
    plexParser(pPlexElt);

    // Check that the plex is connected.
    if(! rParsedPlex.plexIsConnected())
      {
	throw parsedPlexNotConnXcpt(pPlexNode);
      }

    // Recognize the plex, but without doing any of the usual
    // initializations of the resulting cptPlexFamily.
    cptPlexFamily* pFamily = 0;
    cpx::plexIso recognitionIsos;
    rClxUnit.recognize.unify(rParsedPlex,
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
    // with the paradigm of the cptPlexFamily.
    rParsedPlex = parserPlex(nameToMolNdx,
			     pFamily->getParadigm());

    return pFamily;
  }

  cptPlexFamily*
  recognizePlexElt(xmlpp::Element* pPlexElt,
		   parserPlex& rParsedPlex,
		   cml::cmlUnit& rCmlUnit,
		   clxUnit& rClxUnit)
  {
    // Parse the plex.
    parserPlex thePlex;
    parsePlex plexParser(rCmlUnit,
			 thePlex);
    plexParser(pPlexElt);

    // Recognize the plex, but without doing any of the usual
    // initializations of the resulting cptPlexFamily.
    cpx::plexIso recognitionIsos;
    cptPlexFamily* pFamily = rClxUnit.recognize(thePlex,
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
      cml::cmlUnit& rCmlUnit;
      cpt::cptUnit& rCptUnit;
      cptPlexQueries* pQuery;
      cml::cptModMol* pMol;
      cpx::molSpec instanceNdx;
    public:
      parseModSiteRefQuery(cml::cmlUnit& refCmlUnit,
			   cpt::cptUnit& refCptUnit,
			   cptPlexQueries* pAndPlexQueries,
			   cml::cptModMol* pModMol,
			   cpx::molSpec molInstanceNdx) :
	rCmlUnit(refCmlUnit),
	rCptUnit(refCptUnit),
	pQuery(pAndPlexQueries),
	pMol(pModMol),
	instanceNdx(molInstanceNdx)
      {}

      void
      operator()(xmlpp::Node* pModSiteRefNode) const throw(utl::xcpt)
      {
	xmlpp::Element* pModSiteRefElt
	  = utl::dom::mustBeElementPtr(pModSiteRefNode);
    
	// Get the modification site name.
	std::string modSiteName
	  = utl::dom::mustGetAttrString(pModSiteRefElt,
					cml::eltName::modSiteRef_nameAttr);

	// Look up the modification site on the mol, getting the
	// index of the modification site.
	int modSiteNdx = pMol->mustGetModSiteNdx(modSiteName,
						 pModSiteRefElt);

	// Get the modification name.
	xmlpp::Element* pModRefElt
	  = utl::dom::mustGetUniqueChild(pModSiteRefNode,
					 cml::eltName::modRef);
	std::string modName
	  = utl::dom::mustGetAttrString(pModRefElt,
					cml::eltName::modRef_nameAttr);

	// Look up the modification.
	const cpx::modification* pMod
	  = rCmlUnit.mustGetMod(modName,
				pModRefElt);
	
	// Construct query of modMol's state and register for memory
	// management.
	cpx::modMolStateQuery* pMolQuery
	  = new cpx::modMolStateQuery(modSiteNdx,
				      pMod);
	rCptUnit.addQuery(pMolQuery);

	// Construct query of modMol instance's state and register for
	// memory management.
	cpx::molStatePlexQuery<cptPlexSpecies, cptOmniPlex>* pPlexQuery
	  = new cpx::molStatePlexQuery<cptPlexSpecies, cptOmniPlex>
	  (instanceNdx,
	   *pMolQuery);
	rCptUnit.addQuery(pPlexQuery);

	// Add the query to the conjunction of queries to be applied to
	// species of complexes appearing in this stream.
	pQuery->addQuery(pPlexQuery);
      }
    };

    // Parses a single mod-mol-instance-ref element, adding the
    // resulting query to a given, overall query for a plex.
    class parseModMolInstanceStateQuery :
      public std::unary_function<xmlpp::Node*, void>
    {
      cptPlexQueries* pQuery;
      const parserPlex& rPlex;
      cml::cmlUnit& rCmlUnit;
      cpt::cptUnit& rCptUnit;
    
    public:
      parseModMolInstanceStateQuery(cptPlexQueries* pAndPlexQueries,
				    const parserPlex& rParsedPlex,
				    cml::cmlUnit& refCmlUnit,
				    cpt::cptUnit& refCptUnit) :
	pQuery(pAndPlexQueries),
	rPlex(rParsedPlex),
	rCmlUnit(refCmlUnit),
	rCptUnit(refCptUnit)
      {}

      void
      operator()(xmlpp::Node* pModMolInstanceRefNode) const
	throw(utl::xcpt)
      {
	xmlpp::Element* pModMolInstanceRefElt
	  = utl::dom::mustBeElementPtr(pModMolInstanceRefNode);

	std::string instanceName
	  = utl::dom::mustGetAttrString
	  (pModMolInstanceRefElt,
	   eltName::modMolInstanceRef_nameAttr);

	// Get the instance index.  Each filter clause needs a cpx::molSpec,
	// which is an instance index together with the index of the binding
	// on the instance mol.
	int instanceNdx = rPlex.getMolNdxByName(instanceName);
	if(0 > instanceNdx)
	  throw unkMolInstXcpt(instanceName,
			       pModMolInstanceRefElt);

	// Make sure the mol is a mod-mol.
	cml::cptModMol* pModMol
	  = cml::mustBeModMol(rPlex.mols[instanceNdx],
			      pModMolInstanceRefElt);

	// Process the modification map, which gives the state
	// of this instance mol.  This is a map from binding sites to
	// shapes.
	xmlpp::Element* pModMapElt
	  = utl::dom::mustGetUniqueChild(pModMolInstanceRefElt,
					 cml::eltName::modMap);
	xmlpp::Node::NodeList modSiteRefs
	  = pModMapElt->get_children(cml::eltName::modSiteRef);
	std::for_each(modSiteRefs.begin(),
		      modSiteRefs.end(),
		      parseModSiteRefQuery(rCmlUnit,
					   rCptUnit,
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
				    cptPlexQueries* pQuery,
				    const parserPlex& rPlex,
				    cml::cmlUnit& rCmlUnit,
				    cpt::cptUnit& rCptUnit)
      throw(utl::xcpt)
    {
      xmlpp::Element* pInstanceStatesElt
	= utl::dom::mustBeElementPtr(pInstanceStatesNode);
	  
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
						  rCmlUnit,
						  rCptUnit));
    }
  }

  // This routine is really just the place where new parsers for
  // different kinds of mol state specifications could be installed.
  void
  parseInstanceStateQueries(xmlpp::Node* pInstanceStatesNode,
			    cptPlexQueries* pQuery,
			    const parserPlex& rParserPlex,
			    cml::cmlUnit& rCmlUnit,
			    cpt::cptUnit& rCptUnit)
    throw(utl::xcpt)
  {
    // Parse mod-mol state queries, adding them to the overall plex query.
    parseModMolInstanceStateQueries(pInstanceStatesNode,
				    pQuery,
				    rParserPlex,
				    rCmlUnit,
				    rCptUnit);

    // This is where parsers for new kinds of mol state queries could
    // be installed.
  }

  cptPlexFamily*
  parsePlexClass::
  operator()(xmlpp::Node* pParentNode) const
    throw(utl::xcpt)
  {
    xmlpp::Element* pParentElt
      = utl::dom::mustBeElementPtr(pParentNode);

    // Parse the complex.
    xmlpp::Element* pPlexElt
      = utl::dom::mustGetUniqueChild(pParentElt,
				     eltName::plex);
    cptPlexFamily* pFamily = recognizePlexElt(pPlexElt,
					   rPlex,
					   rCmlUnit,
					   rClxUnit);

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
				  rCmlUnit,
				  rCptUnit);
      }

    return pFamily;
  }

  namespace
  {
    // Parses a single mod-site-ref/modification pair, converting it
    // into an entry in a map from modification site name to modification.
    class parseModSiteRef :
      public std::unary_function
    <xmlpp::Node*, std::pair<std::string, const cpx::modification*> >
    {
      cml::cmlUnit& rCmlUnit;
    
    public:
      parseModSiteRef(cml::cmlUnit& refCmlUnit) :
	rCmlUnit(refCmlUnit)
      {}

      std::pair<std::string, const cpx::modification*>
      operator()(xmlpp::Node* pModSiteRefNode) const
	throw(utl::xcpt)
      {
	xmlpp::Element* pModSiteRefElt
	  = utl::dom::mustBeElementPtr(pModSiteRefNode);

	// Get the modification site name.
	std::string modSiteName
	  = utl::dom::mustGetAttrString(pModSiteRefElt,
					cml::eltName::modSiteRef_nameAttr);

	// Get the modification element, and from it, the name of the
	// modification.
	xmlpp::Element* pModRefElt
	  = utl::dom::mustGetUniqueChild(pModSiteRefElt,
					 cml::eltName::modRef);
	std::string modName
	  = utl::dom::mustGetAttrString(pModRefElt,
					cml::eltName::modRef_nameAttr);

	// Look up the modification in cml::cmlUnit.
	return std::make_pair(modSiteName,
			      rCmlUnit.mustGetMod(modName));
      }
    };

    class parseModMolInstanceState :
      public std::unary_function<xmlpp::Node*, void>
    {
      std::vector<cpx::molParam>& rParams;
      const parserPlex& rPlex;
      const cptPlexFamily* pFamily;
      cml::cmlUnit& rCmlUnit;

    public:
      parseModMolInstanceState(std::vector<cpx::molParam>& rMolParams,
			       const parserPlex& rParserPlex,
			       const cptPlexFamily* pPlexFamily,
			       cml::cmlUnit& refCmlUnit) :
	rParams(rMolParams),
	rPlex(rParserPlex),
	pFamily(pPlexFamily),
	rCmlUnit(refCmlUnit)
      {}

      void
      operator()(xmlpp::Node* pModMolInstanceRefNode) const
	throw(utl::xcpt)
      {
	// Make sure the mod-mol-instance-ref is an element; probably
	// unnecessarily dynamic cast.
	xmlpp::Element* pModMolInstanceRefElt
	  = utl::dom::mustBeElementPtr(pModMolInstanceRefNode);

	// Get the instance name of the mod-mol instance.
	std::string instanceName
	  = utl::dom::mustGetAttrString(pModMolInstanceRefElt,
					eltName::modMolInstanceRef_nameAttr);

	// Get the index of the mod-mol instance.
	int instanceNdx
	  = rPlex.getMolNdxByName(instanceName);

	// Get the mod-mol itself and verify that it is a mod-mol.
	const cptPlex& rParadigm = pFamily->getParadigm();
	cml::cptModMol* pModMol
	  = cml::mustBeModMol(rParadigm.mols[instanceNdx],
			      pModMolInstanceRefElt);

	// Target for parsing mod-map elements.  The mod mol
	// can intern this to to generate the corresponding modMolState*.
	std::map<std::string, const cpx::modification*> modMap;

	// Parse the modifcation map.
	xmlpp::Element* pModMapElt
	  = utl::dom::mustGetUniqueChild(pModMolInstanceRefElt,
					 cml::eltName::modMap);
	xmlpp::Node::NodeList modSiteRefs
	  = pModMapElt->get_children(cml::eltName::modSiteRef);
	std::transform(modSiteRefs.begin(),
		       modSiteRefs.end(),
		       std::inserter(modMap,
				     modMap.begin()),
		       parseModSiteRef(rCmlUnit));

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
			      std::vector<cpx::molParam>& rMolParams,
			      const parserPlex& rParserPlex,
			      const cptPlexFamily* pPlexFamily,
			      cml::cmlUnit& rCmlUnit)
    {
      xmlpp::Element* pInstanceStatesElt
	= utl::dom::mustBeElementPtr(pInstanceStatesNode);
	  
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
					     rCmlUnit));
    }
  }

  void
  parseInstanceStates(xmlpp::Node* pInstanceStatesNode,
		      std::vector<cpx::molParam>& rMolParams,
		      const parserPlex& rParserPlex,
		      const cptPlexFamily* pPlexFamily,
		      cml::cmlUnit& rCmlUnit)
    throw(utl::xcpt)
  {
    // Parse mod-mol states, replacing the states in rMolParams.
    parseModMolInstanceStates(pInstanceStatesNode,
			      rMolParams,
			      rParserPlex,
			      pPlexFamily,
			      rCmlUnit);

    // This is where new parsers for new kinds of mol states for new kinds
    // of mols could be installed.
  }

  cptPlexSpecies*
  parsePlexSpecies(xmlpp::Element* pParentElt,
		   cpt::cptUnit& rCptUnit,
		   cml::cmlUnit& rCmlUnit,
		   clxUnit& rClxUnit)
    throw(utl::xcpt)
  {
    // Get the cptPlexFamily and arrange for remapping mol instance names
    // to indices in the paradigm plex.
    xmlpp::Element* pPlexElt
      = utl::dom::mustGetUniqueChild(pParentElt,
				     eltName::plex);
    parserPlex thePlex;
    cptPlexFamily* pPlexFamily
      = recognizePlexElt(pPlexElt,
			 thePlex,
			 rCmlUnit,
			 rClxUnit);

    // If no instance states are specified, then we use the default molParams.
    std::vector<cpx::molParam> molParams
      = pPlexFamily->makeDefaultMolParams();

    // Parse optional instanceStates of the instance mols. These will replace
    // some of the molParams.
    xmlpp::Element* pInstanceStatesElt
      = utl::dom::getOptionalChild(pParentElt,
				   eltName::instanceStates);
    if(pInstanceStatesElt)
      {
	parseInstanceStates(pInstanceStatesElt,
			    molParams,
			    thePlex,
			    pPlexFamily,
			    rCmlUnit);
      }
	  
    // Construct the new cptPlexSpecies from the molParams using the
    // parsed plex's family.
    return pPlexFamily->getMember(molParams);
  }

  cptPlexSpecies*
  parseExplicitPlexSpecies::
  operator()(xmlpp::Node* pPlexSpeciesNode) const
    throw(utl::xcpt)
  {
    xmlpp::Element* pPlexSpeciesElt
      = utl::dom::mustBeElementPtr(pPlexSpeciesNode);

    // Get the name of the plex species.
    std::string plexSpeciesName
      = utl::dom::mustGetAttrString(pPlexSpeciesElt,
				    eltName::plexSpecies_nameAttr);

    // Parse the plexSpecies.
    cptPlexSpecies* pPlexSpecies
      = parsePlexSpecies(pPlexSpeciesElt,
			 rCptUnit,
			 rCmlUnit,
			 rClxUnit);

    // Insert the plexSpecies into cptUnit's catalog of named species.
    rCptUnit.mustNameSpecies(plexSpeciesName,
			     pPlexSpecies,
			     pPlexSpeciesNode);

    return pPlexSpecies;
  }

  cptPlexSpecies*
  parseTaggedPlexSpecies::
  operator()(xmlpp::Node* pTaggedPlexSpeciesNode) const
    throw(utl::xcpt)
  {
    xmlpp::Element* pTaggedPlexSpeciesElt
      = utl::dom::mustBeElementPtr(pTaggedPlexSpeciesNode);

    // Parse the plexSpecies.
    cptPlexSpecies* pPlexSpecies
      = parsePlexSpecies(pTaggedPlexSpeciesElt,
			 rCptUnit,
			 rCmlUnit,
			 rClxUnit);

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
      public std::unary_function<xmlpp::Node*, cpx::siteToShapeMap::value_type>
    {
      const parserPlex& rPlex;

    public:
      parseMolInstanceRef(const parserPlex& rParserPlex) :
	rPlex(rParserPlex)
      {}

      result_type
      operator()(argument_type pMolInstanceRefNode) const
	throw(utl::xcpt)
      {
	// Check that the node is an element node, possibly unnecessarily.
	xmlpp::Element* pMolInstanceRefElt
	  = utl::dom::mustBeElementPtr(pMolInstanceRefNode);
    
	// Get the mol instance name.
	std::string molInstanceName
	  = utl::dom::mustGetAttrString
	  (pMolInstanceRefElt,
	   eltName::molInstanceRef_nameAttr);

	// Get the corresponding mol index, part of the cpx::siteSpec, in order
	// to make an entry in the allostery map of the plex definition.
	int molNdx
	  = rPlex.mustGetMolNdxByName(pMolInstanceRefElt,
				      molInstanceName);

	// Get the corresponding mol, in order to look up binding site by
	// name.
	cml::cptMol* pMol = rPlex.mols[molNdx];

	// Get the child element that gives the site name.
	xmlpp::Element* pBindingSiteRefElt
	  = utl::dom::mustGetUniqueChild(pMolInstanceRefElt,
					 cml::eltName::bindingSiteRef);

	// Get the name of the binding site.
	std::string bindingSiteName
	  = utl::dom::mustGetAttrString
	  (pBindingSiteRefElt,
	   cml::eltName::bindingSiteRef_nameAttr);

	// Get the index of the binding site, part of the cpx::siteSpec.
	int bindingSiteNdx
	  = pMol->mustFindSite(bindingSiteName,
			       pBindingSiteRefElt);

	// Get the actual binding site, to look up the shape.
	cml::cptMol& rMol = *pMol;
	cml::cptBndSite& rBindingSite
	  = rMol[bindingSiteNdx];
    
	// Get the child element that gives the shape of the binding
	// site when in this subcomplex.
	xmlpp::Element* pSiteShapeRefElt
	  = utl::dom::mustGetUniqueChild(pBindingSiteRefElt,
					 cml::eltName::siteShapeRef);

	// Get the name of the site shape.
	std::string siteShapeName
	  = utl::dom::mustGetAttrString(pSiteShapeRefElt,
					cml::eltName::siteShapeRef_nameAttr);

	// Look up the site shape.
	const cpx::siteShape* pSiteShape
	  = rBindingSite.mustGetShape(pMol,
				      siteShapeName,
				      pSiteShapeRefElt);

	// Return entry for allostery map.
	return std::make_pair(cpx::siteSpec(molNdx,
					   bindingSiteNdx),
			      pSiteShape);
      }
    };
  }

  void
  parseAllostericSites::
  operator()(xmlpp::Element* pAlloSitesElt) const
    throw(utl::xcpt)
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
    throw(utl::xcpt)
  {
    // Parse the plex class.
    parserPlex parsedPlex;

    cptPlexQueries* pAndPlexQueries
      = new cptPlexQueries();
    rCptUnit.addQuery(pAndPlexQueries);

    parsePlexClass classParser(rCptUnit,
			       rCmlUnit,
			       rClxUnit,
			       parsedPlex,
			       pAndPlexQueries);

    cptPlexFamily* pFamily
      = classParser(pAlloPlexNode);

    // Parse the allosteric sites.
    cpx::siteToShapeMap alloSiteMap;
    parseAllostericSites siteParser(parsedPlex,
				    alloSiteMap);

    // Add the query and alloSiteMap to the cptPlexFamily's alloStateList, which
    // is used to test new plexSpecies as they appear, applying allosteric
    // modifications to those that answer the query.
    pFamily->addAlloQueryAndMap(pAndPlexQueries,
				alloSiteMap);
  }

  void
  parsePlexSpeciesStream::
  operator()(xmlpp::Node* pPlexSpeciesStreamNode) const
    throw(utl::xcpt)
  {
    // Make sure that the node is an element, possibly unnecessarily.
    xmlpp::Element* pPlexSpeciesStreamElt
      = utl::dom::mustBeElementPtr(pPlexSpeciesStreamNode);

    // Get the name of the species stream.
    std::string streamName
      = utl::dom::mustGetAttrString
      (pPlexSpeciesStreamElt,
       eltName::plexSpeciesStream_nameAttr);

    // Parse the plex class, to get a plex, cptPlexFamily, and a state query.
    parserPlex parsedPlex;
    cptPlexQueries* pAndPlexQueries = new cptPlexQueries();
    rCptUnit.addQuery(pAndPlexQueries);
    
    // Is this ever used in an STL algorithm? Maybe just a function?  Why
    // do I bother?
    parsePlexClass classParser(rCptUnit,
			       rCmlUnit,
			       rClxUnit,
			       parsedPlex,
			       pAndPlexQueries);
    cptPlexFamily* pFamily = classParser(pPlexSpeciesStreamElt);

    // Create the dumpable.
    cpt::queryGlobalSpeciesDumpable<cptPlexSpecies>* pDumpable
      = new cpt::queryGlobalSpeciesDumpable<cptPlexSpecies>(streamName,
							    rCptUnit.getCompartmentGraph(),
							    *pAndPlexQueries);

    // Add the dumpable to the cptPlexFamily.
    pFamily->addDumpable(pDumpable);

    // Add the dumpable to moleculizer for deletion and access.
    rCptUnit.addSpeciesDumpable(pDumpable);
  }
}
