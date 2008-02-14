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

#include "clx/parseOmniPlex.hh"
#include "clx/cptPlexFamily.hh"
#include "clx/parsePlex.hh"
#include "clx/cptPlexQueries.hh"
#include "clx/cptOmniStructureQueries.hh"

namespace clx
{
  class parseUnboundSiteQuery :
    public std::unary_function<xmlpp::Node*, void>
  {
    cptOmniStructureQueries* pQueries;
    const parserPlex& rParserPlex;
    clxUnit& rClxUnit;
    
  public:
    parseUnboundSiteQuery(cptOmniStructureQueries* pStructureQueries,
			  const parserPlex& rParsedPlex,
			  clxUnit& refClxUnit) :
      pQueries(pStructureQueries),
      rParserPlex(rParsedPlex),
      rClxUnit(refClxUnit)
    {}

    void
    operator()(xmlpp::Node* pInstanceRefNode) const
      throw(utl::xcpt)
    {
      xmlpp::Element* pInstanceRefElt
	= utl::dom::mustBeElementPtr(pInstanceRefNode);

      // Parse the instance name.
      std::string instanceName
	= utl::dom::mustGetAttrString(pInstanceRefElt,
				      eltName::instanceRef_nameAttr);

      // Convert the instance name to a mol index.
      int molNdx = rParserPlex.mustGetMolNdxByName(pInstanceRefElt,
						   instanceName);

      xmlpp::Element* pSiteRefElt
	= utl::dom::mustGetUniqueChild(pInstanceRefElt,
				       eltName::siteRef);

      // Parse the name of the site that is supposed to be free.
      std::string siteName
	= utl::dom::mustGetAttrString(pSiteRefElt,
				      eltName::siteRef_nameAttr);

      // Ask the mol to convert the site name into a site index.
      const cml::cptMol* pMol = rParserPlex.mols[molNdx];
      int siteNdx = pMol->mustFindSite(siteName,
				       pSiteRefElt);

      // Construct the query, and add it to the clxUnit for memory
      // management.
      cptOmniFreeSiteQuery* pFreeSiteQuery
	= new cptOmniFreeSiteQuery(cpx::siteSpec(molNdx,
						 siteNdx));
      rClxUnit.addStructureQuery(pFreeSiteQuery);

      // Add the free site query to the overall structural query.
      pQueries->addQuery(pFreeSiteQuery);
    }
  };
    
  void
  parseOmniPlex::
  operator()(xmlpp::Node* pParentNode) const
    throw(utl::xcpt)
  {
    // Unify the plex; i.e. find its plexFamily in the database, or create it,
    // but don't initialize the plexFamily (connectToFeatures).  plexFamilies
    // can't be connected to their features until after all omniPlexes have
    // been parsed in this way.
    xmlpp::Element* pPlexElt
      = utl::dom::mustGetUniqueChild(pParentNode,
				     eltName::plex);
    parserPlex parsedPlex;
    cptPlexFamily* pFamily
      = unifyPlexNode(pPlexElt,
		      rCmlUnit,
		      rClxUnit,
		      parsedPlex);

    // Parse the instance states, getting a query.
    cptPlexQueries* pAndPlexQueries
      = new cptPlexQueries();
    
    xmlpp::Element* pInstanceStatesElt
      = utl::dom::getOptionalChild(pParentNode,
				   eltName::instanceStates);
    if(pInstanceStatesElt)
      {
	parseInstanceStateQueries(pInstanceStatesElt,
				  pAndPlexQueries,
				  parsedPlex,
				  rCmlUnit,
				  rCptUnit);
      }

    // Parse optional structure queries.
    //
    // The only structure query for the time being is a test if a particular
    // (free) site on the omni is free in the plex where the omni is found.
    cptOmniStructureQueries* pStructureQueries
      = new cptOmniStructureQueries();
    rClxUnit.addStructureQuery(pStructureQueries);

    xmlpp::Element* pUnboundSitesElt
      = utl::dom::getOptionalChild(pParentNode,
				   eltName::unboundSites);
    if(pUnboundSitesElt)
      {
	xmlpp::Node::NodeList instanceRefNodes
	  = pUnboundSitesElt->get_children(eltName::instanceRef);

	std::for_each(instanceRefNodes.begin(),
		      instanceRefNodes.end(),
		      parseUnboundSiteQuery(pStructureQueries,
					    parsedPlex,
					    rClxUnit));
      }

    // Construct the cptOmniPlex (which also adds it to its plexFamily.)
    //
    // This is so that new plexFamily's can run down the list
    // of omniPlexes in each structural class that the recognizer
    // found in the new plexFamily.  The new plexFamily checks the
    // structural query of each of these omniPlexes, connecting itself
    // to those whose tests it passes.
    cptOmniPlex* pOmni
      = new cptOmniPlex(pFamily,
			pStructureQueries,
			pAndPlexQueries);

    // Register the family has having omniplexes.
    // 
    // This is so recognizer will check for its presence in
    // new plexes.
    rClxUnit.addOmniPlex(pOmni,
			  pParentNode);
  }

  // This could also return the plexFamily with no additional work.
  cptOmniPlex*
  findOmni(xmlpp::Node* pParentNode,
	   cml::cmlUnit& rCmlUnit,
	   clxUnit& rClxUnit,
	   parserPlex& rParsedPlex)
    throw(utl::xcpt)
  {
    // Get the plexFamily of the omniPlex under pParentNode.
    xmlpp::Element* pPlexElt
      = utl::dom::mustGetUniqueChild(pParentNode,
				     eltName::plex);

    // Here, the cptPlexFamily is returned, but discarded.
    unifyPlexNode(pPlexElt,
		  rCmlUnit,
		  rClxUnit,
		  rParsedPlex);

    // With the reorganization of omniPlexes, this lookup is now separate
    // from all the above.  Previously, the plexFamily we just found did this.
    return rClxUnit.mustGetOmniForNode(pParentNode);
  }

  cptOmniPlex*
  parseAllostericOmni::
  operator()(xmlpp::Node* pParentNode) const
    throw(utl::xcpt)
  {
    // Find omniPlex parsed earlier.
    parserPlex parsedPlex;
    cptOmniPlex* pOmni
      = findOmni(pParentNode,
		 rCmlUnit,
		 rClxUnit,
		 parsedPlex);

    // Parse allosteric sites.
    //
    // The allosteric modifications are installed in the omni's
    // siteToShapeMap.
    parseAllostericSites alloSitesParser(parsedPlex,
					 pOmni->getSiteToShapeMap());
    xmlpp::Element* pAlloSitesElt
      = utl::dom::mustGetUniqueChild(pParentNode,
				     eltName::allostericSites);
    alloSitesParser(pAlloSitesElt);

    return pOmni;
  }

  void
  parseOmniSpeciesStream::
  operator()(xmlpp::Node* pOmniSpeciesStreamNode) const
    throw(utl::xcpt)
  {
    xmlpp::Element* pOmniSpeciesStreamElt
      = utl::dom::mustBeElementPtr(pOmniSpeciesStreamNode);

    // Get the name of the species stream.
    std::string streamName
      = utl::dom::mustGetAttrString
      (pOmniSpeciesStreamElt,
       eltName::omniSpeciesStream_nameAttr);

    // Find omniPlex parsed earlier.
    parserPlex parsedPlex;
    cptOmniPlex* pOmni
      = findOmni(pOmniSpeciesStreamNode,
		 rCmlUnit,
		 rClxUnit,
		 parsedPlex);

    // Construct species dumpable to attach to omniplex.
    // Add it to cptUnit for memory management and lookup.
    cpt::multiGlobalSpeciesDumpable<cptPlexSpecies>* pDumpable
      = new cpt::multiGlobalSpeciesDumpable<cptPlexSpecies>(streamName,
							    rCptUnit.getCompartmentGraph());
    rCptUnit.addSpeciesDumpable(pDumpable);

    // Attach dumpable to omniplex, where it will be told of all
    // new species satisfying the omniplex's queries.
    pOmni->getSubPlexFeature()->setDumpable(pDumpable);
  }
}
