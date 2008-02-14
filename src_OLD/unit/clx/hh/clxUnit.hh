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

#ifndef CLX_CLXUNIT_H
#define CLX_CLXUNIT_H

#include <algorithm>
#include "utl/platform.hh"
#include "cpx/plexQuery.hh"
#include "cpx/omniStructureQuery.hh"
#include "cpx/omniPlex.hh"
#include "cpt/cptUnit.hh"
#include "cpt/cptEltName.hh"
#include "cpt/cptApp.hh"
#include "cpt/unitsMgr.hh"
#include "cml/cptMol.hh"
#include "cml/cmlUnit.hh"
#include "clx/clxEltName.hh"
#include "clx/cptPlexSpecies.hh"
#include "clx/cptRecognizer.hh"
#include "clx/cptPlexFamily.hh"
#include "clx/cptOmniPlex.hh"

namespace clx
{
  class clxUnit : 
    public cpt::unit
  {
  public:
    typedef cpx::plexQuery<cptPlexSpecies,
			   cptOmniPlex> plexQueryType;
  
  private:
    friend class cptRecognizer;
    
    cpt::cptUnit& rCptUnit;
    cml::cmlUnit& rCmlUnit;

    cpx::knownBindings<cml::cptMol, fnd::feature<cpx::cxBinding<cptPlexSpecies, cptPlexFamily> > > bindingFeatures;

    // Another unit that uses plexUnit facilities can register an Xpath, so
    // that the plexUnit can find the other unit's omniplexes.
    std::vector<std::string> omniXpaths;

    // The set of plexFamilies that have associated omniPlexes.
    // This set is for use in the first step of determining
    // which omniPlexes are found in a new plexFamily (i.e. structural
    // class of complexes).  The omniPlexes associated to these
    // plexFamilies may have structural queries (for now, free site queries)
    // that must also be satisfied by the new structural class (plexFamily)
    // in order for the plexFamily to "display" the omniPlex's feature
    // or to be influenced by allosteric properties given by the omniPlex.
    //
    // This is an std::set because several omniPlex specification may
    // point to the same plexFamily, and we want only one entry here
    // per plexFamily to minimize the number of structural comparisons
    // that must be made.
    //
    // Right now, it looks like every plexFamily has a reference to this set.
    std::set<cptPlexFamily*> omniPlexFamilies;

    // Map from omniplexes (frequently parsed by the plexUnit for its client
    // units) to the nodes that were parsed to produce them (by the plexUnit
    // for its client units.  These are looked up later for the clients
    // using this map.
    std::map<xmlpp::Node*,
	     cptOmniPlex*> nodeToOmni;

    // Memory management for plexQueries.
    //
    // PlexQueries may refer to one another: compound queries use the queries
    // that they combine.  Value-based treatment doesn't seem possible, since
    // the essential character of these things is virtual functions.  Hence,
    // we intern all plexQueries here (they are all derived from input
    // elements) and collect them all at the end of time.
    utl::autoVector<plexQueryType> plexQueries;

    // Memory management for omniStructureQueries.
    //
    // Comments above about state queries (plexQueries) also apply here.
    utl::autoVector<fnd::query<cpx::omniStructureQueryArg<cptPlex> > > structureQueries;

  public:
    clxUnit(cpt::cptApp& rCptApp,
	    cpt::cptUnit& refCptUnit,
	    cml::cmlUnit& refCmlUnit);

    /*! \name Database of binding features.

    Each binding feature is connected to a %pair of "structural sites."
    A "structural site" is a a %pair of a mol and a site index on that
    mol.

    Any plexFamily in which the two given structural sites are bound
    together retains a pointer to the bindingFeature associated to the
    %pair of structural sites.  The plexFamily uses the bindingFeature
    to inform the family of decompositions of the binding when a new
    member species is created in the plexFamily.  The decomposition
    family's reaction generator can then create the decomposition
    reactions involving the new species.  */
    //@{
    /*! \brief Add binding feature to database.

    Used in the dimer command, where decomposition reaction families
    are created. The decompFam for bindings of the given "structural
    sites" is notified through the bindingFeature when a new species
    of complex appears that contains the given binding.*/
    fnd::feature<cpx::cxBinding<cptPlexSpecies, cptPlexFamily> >*
    addBindingFeature(cml::cptMol* pLeftMol,
		      int leftMolSiteSpec,
		      cml::cptMol* pRightMol,
		      int rightMolSiteSpec);

    /*! \brief Recover binding feature for %pair of structural sites.

    This is used to connect a plexFamily to its binding features.
    Doing so enables the plexFamily to notify the reaction families
    that are interested in bindings (only decompFam ?) when a new
    member species has been added to the plexFamily. */
    fnd::feature<cpx::cxBinding<cptPlexSpecies, cptPlexFamily> >*
    findBindingFeature(cml::cptMol* pLeftMol,
		       int leftMolSiteSpec,
		       cml::cptMol* pRightMol,
		       int rightMolSiteSpec);


    // When an omniplex is parsed (possibly at the behest of another unit)
    // it should be added here.
    void
    addOmniPlex(cptOmniPlex* pOmniPlex,
		xmlpp::Node* pParentNode)
      throw(utl::xcpt);

    // For use by plexFamily::connectToFeatures.  This allows one to determine
    // which omniplexes (i.e. distinguished subcomplexes) appear in a new
    // structural family of complex species.
    const std::set<cptPlexFamily*>&
    getOmniPlexFamilies(void) const
    {
      return omniPlexFamilies;
    }

    // The plexUnit has to parse omniPlexes for other units, so that
    // they are known before plex species, etc. are parsed.  This is
    // how the plexUnit makes the parsing result available to those
    // other units.
    //
    // Looks up the omniPlex parsed for a particular node (there is no
    // "omniPlex node," but an omniplex's nodes are all beneath the same
    // parent, which is the node argument here.)  Returns 0 if none is found,
    // which should be an internal error.
    cptOmniPlex*
    getOmniForNode(xmlpp::Node* pParentNode) const;

    // Looks up the omniPlex parsed for a particular node (there is no
    // "omniPlex node," but an omniplex's nodes are all beneath the same
    // parent, which is the node argument here.)  Throws an "internal
    // exception" if no omniplex is found.
    cptOmniPlex*
    mustGetOmniForNode(xmlpp::Node* pParentNode) const
      throw(utl::xcpt);

    cptRecognizer recognize;

    void
    addPlexQuery(plexQueryType* pQuery)
    {
      plexQueries.push_back(pQuery);
    }

    void
    addStructureQuery(fnd::query<cpx::omniStructureQueryArg<cptPlex> >* pQuery)
    {
      structureQueries.push_back(pQuery);
    }

    void
    addOmniXpath(const std::string& rXpath)
    {
      omniXpaths.push_back(rXpath);
    }

    void
    parseDomInput(xmlpp::Element* pRootElt,
		  xmlpp::Element* pModelElt,
		  xmlpp::Element* pStreamsElt,
		  xmlpp::Element* pEventsElt)
      throw(utl::xcpt);
  
    void
    insertStateElts(xmlpp::Element* pRootElt)
      throw(std::exception);

    // Here, all the dumpables and reactions are notified of all
    // the existing plex species all at once.
    void
    prepareToRun(xmlpp::Element* pRootElt,
		 xmlpp::Element* pModelElt,
		 xmlpp::Element* pStreamsElt,
		 xmlpp::Element* pEventsElt)
      throw(utl::xcpt);

    void
    prepareToDump(xmlpp::Element* pRootElt,
		  xmlpp::Element* pModelElt,
		  xmlpp::Element* pStreamsElt,
		  xmlpp::Element* pEventsElt,
		  xmlpp::Element* pTaggedSpeciesElement)
      throw(utl::xcpt);

    void
    prepareToContinue(xmlpp::Element* pRootElt,
		      xmlpp::Element* pModelElt,
		      xmlpp::Element* pStreamsElt,
		      xmlpp::Element* pEventsElt,
		      std::map<std::string, std::string>& rTagToName,
		      xmlpp::Element* pTaggedSpeciesElement)
      throw(utl::xcpt);
  };
}

#endif // CLX_CLXUNIT_H
