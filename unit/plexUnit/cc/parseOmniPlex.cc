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

#include "mzr/dumpable.hh"
#include "plex/parseOmniPlex.hh"
#include "plex/plexFamily.hh"
#include "plex/parsePlex.hh"

namespace plx
{
  class parseUnboundSiteQuery :
    public std::unary_function<xmlpp::Node*, void>
  {
    andOmniStructureQueries* pQueries;
    const parserPlex& rParserPlex;
    plexUnit& rPlexUnit;
    
  public:
    parseUnboundSiteQuery(andOmniStructureQueries* pStructureQueries,
			  const parserPlex& rParsedPlex,
			  plexUnit& refPlexUnit) :
      pQueries(pStructureQueries),
      rParserPlex(rParsedPlex),
      rPlexUnit(refPlexUnit)
    {}

    void
    operator()(xmlpp::Node* pInstanceRefNode) const
      throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pInstanceRefElt
	= domUtils::mustBeElementPtr(pInstanceRefNode);

      // Parse the instance name.
      std::string instanceName
	= domUtils::mustGetAttrString(pInstanceRefElt,
				      eltName::instanceRef_nameAttr);

      // Convert the instance name to a mol index.
      int molNdx = rParserPlex.mustGetMolNdxByName(pInstanceRefElt,
						   instanceName);

      xmlpp::Element* pSiteRefElt
	= domUtils::mustGetUniqueChild(pInstanceRefElt,
				       eltName::siteRef);

      // Parse the name of the site that is supposed to be free.
      std::string siteName
	= domUtils::mustGetAttrString(pSiteRefElt,
				      eltName::siteRef_nameAttr);

      // Ask the mol to convert the site name into a site index.
      const bnd::mol* pMol = rParserPlex.mols[molNdx];
      int siteNdx = pMol->mustGetSiteNdxByName(pSiteRefElt,
					       siteName);

      // Construct the query, and add it to the plexUnit for memory
      // management.
      omniFreeSiteQuery* pFreeSiteQuery
	= new omniFreeSiteQuery(plexSiteSpec(molNdx,
					     siteNdx));
      rPlexUnit.addStructureQuery(pFreeSiteQuery);

      // Add the free site query to the overall structural query.
      pQueries->addQuery(pFreeSiteQuery);
    }
  };
    
  void
  parseOmniPlex::
  operator()(xmlpp::Node* pParentNode) const
    throw(mzr::mzrXcpt)
  {
    // Unify the plex; i.e. find its plexFamily in the database, or create it,
    // but don't initialize the plexFamily (connectToFeatures).  plexFamilies
    // can't be connected to their features until after all omniPlexes have
    // been parsed in this way.
    xmlpp::Element* pPlexElt
      = domUtils::mustGetUniqueChild(pParentNode,
				     eltName::plex);
    parserPlex parsedPlex;
    plexFamily* pFamily
      = unifyPlexNode(pPlexElt,
		      rMolUnit,
		      rPlexUnit,
		      parsedPlex);

    // Parse the instance states, getting a query.
    andPlexQueries* pAndPlexQueries
      = new andPlexQueries(rPlexUnit);
    
    xmlpp::Element* pInstanceStatesElt
      = domUtils::getOptionalChild(pParentNode,
				   eltName::instanceStates);
    if(pInstanceStatesElt)
      {
	parseInstanceStateQueries(pInstanceStatesElt,
				  pAndPlexQueries,
				  parsedPlex,
				  rMolUnit,
				  rPlexUnit);
      }

    // Parse optional structure queries.
    //
    // The only structure query for the time being is a test if a particular
    // (free) site on the omni is free in the plex where the omni is found.
    andOmniStructureQueries* pStructureQueries
      = new andOmniStructureQueries();
    rPlexUnit.addStructureQuery(pStructureQueries);

    xmlpp::Element* pUnboundSitesElt
      = domUtils::getOptionalChild(pParentNode,
				   eltName::unboundSites);
    if(pUnboundSitesElt)
      {
	xmlpp::Node::NodeList instanceRefNodes
	  = pUnboundSitesElt->get_children(eltName::instanceRef);

	std::for_each(instanceRefNodes.begin(),
		      instanceRefNodes.end(),
		      parseUnboundSiteQuery(pStructureQueries,
					    parsedPlex,
					    rPlexUnit));
      }

    // Construct the omniPlex; this constructor also adds the omniPlex
    // to its plexFamily.
    //
    // This is so that new plexFamily's can run down the list
    // of omniPlexes in each structural class that the recognizer
    // found in the new plexFamily.  The new plexFamily checks the
    // structural query of each of these omniPlexes, connecting itself
    // to those whose tests it passes.
    new omniPlex(pParentNode,
		 pFamily,
		 pStructureQueries,
		 pAndPlexQueries);

    // Register the family has having omniplexes.
    // 
    // This is so recognizer will check for its presence in
    // new plexes.
    rPlexUnit.addOmniPlexFamily(pFamily);
  }

  omniPlex*
  findOmni(xmlpp::Node* pParentNode,
	   bnd::molUnit& rMolUnit,
	   plexUnit& rPlexUnit,
	   parserPlex& rParsedPlex)
    throw(mzr::mzrXcpt)
  {
    // Get the plexFamily of the omniPlex under pParentNode.
    xmlpp::Element* pPlexElt
      = domUtils::mustGetUniqueChild(pParentNode,
				     eltName::plex);
    plexFamily* pFamily
      = unifyPlexNode(pPlexElt,
		      rMolUnit,
		      rPlexUnit,
		      rParsedPlex);

    // Run through the plexFamily's omniPlexes, searching for
    // the one associated to pParentNode.
    return pFamily->mustFindOmniForNode(pParentNode);
  }

  omniPlex*
  parseAllostericOmni::
  operator()(xmlpp::Node* pParentNode) const
    throw(mzr::mzrXcpt)
  {
    // Find omniPlex parsed earlier.
    parserPlex parsedPlex;
    omniPlex* pOmni
      = findOmni(pParentNode,
		 rMolUnit,
		 rPlexUnit,
		 parsedPlex);

    // Parse allosteric sites.
    //
    // The allosteric modifications are installed in the omni's
    // siteToShapeMap.
    parseAllostericSites alloSitesParser(parsedPlex,
					 pOmni->getSiteToShapeMap());
    xmlpp::Element* pAlloSitesElt
      = domUtils::mustGetUniqueChild(pParentNode,
				     eltName::allostericSites);
    alloSitesParser(pAlloSitesElt);

    return pOmni;
  }

  void
  parseOmniSpeciesStream::
  operator()(xmlpp::Node* pOmniSpeciesStreamNode) const
    throw(mzr::mzrXcpt)
  {
    xmlpp::Element* pOmniSpeciesStreamElt
      = domUtils::mustBeElementPtr(pOmniSpeciesStreamNode);

    // Get the name of the species stream.
    std::string streamName
      = domUtils::mustGetAttrString
      (pOmniSpeciesStreamElt,
       eltName::omniSpeciesStream_nameAttr);

    // Find omniPlex parsed earlier.
    parserPlex parsedPlex;
    omniPlex* pOmni
      = findOmni(pOmniSpeciesStreamNode,
		 rMolUnit,
		 rPlexUnit,
		 parsedPlex);

    // Construct species dumpable to attach to omniplex.
    // Add it to mzrUnit for memory management and lookup.
    mzr::multiSpeciesDumpable* pDumpable
      = new mzr::multiSpeciesDumpable(streamName);
    rMzrUnit.addSpeciesDumpable(pDumpable);

    // Attach dumpable to omniplex, where it will be told of all
    // new species satisfying the omniplex's queries.
    pOmni->getSubPlexFeature()->setDumpable(pDumpable);
  }
}
