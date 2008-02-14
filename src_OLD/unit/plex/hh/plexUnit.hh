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

#ifndef PLEXUNIT_H
#define PLEXUNIT_H

#include <algorithm>
#include "utl/platform.hh"
#include "cpx/plexQuery.hh"
#include "cpx/omniStructureQuery.hh"
#include "cpx/omniPlex.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/moleculizer.hh"
#include "mzr/unitsMgr.hh"
#include "mol/mzrMol.hh"
#include "mol/molUnit.hh"
#include "plex/plexEltName.hh"
#include "plex/mzrPlexSpecies.hh"
#include "plex/mzrRecognizer.hh"
#include "plex/mzrPlexFamily.hh"
#include "plex/mzrOmniPlex.hh"

namespace plx
{
  /*! \defgroup plexGroup The plex unit.
    \ingroup unitsGroup
    \brief Provides complexes and subcomplexes.

    This unit provides plexSpecies, which are species of protein
    complexes, and plexFamily, the single parameterized species at this
    point in development.  A plexFamily is a family of complexes that
    all have the same structure; they are differentiated from one
    another by the internal states of their mols.  A plexFamily is a
    highly embellished STL map from (vectors of) mol states to precise
    species of complex.

    This unit also provides the recognizer.  This mechanism determines
    the structural class of a complex constructed from scratch by a
    reaction.  This is a fundamental part of the automatic manipulation
    of new species of complexes.  Once the structure is determined,
    the plexFamily for the complex is known.  The plexFamily is based
    on an STL map from mol states to precise species of complex, and
    a map lookup in this map is the final step in recognizing a species
    of complex.

    This unit gives the underpinning for selecting plex species for
    dumping complexes or not based on the state of their constituent
    mols.  Basic mols don't have any state, so that these underpinnings
    don't get useful until you load a unit that defines an interesting
    class of mols (e.g. phos) and queries to determine their interesting
    states.  */

  /*! \file plexUnit.hh
    \ingroup plexGroup
    \brief Defines plexUnit, a unit providing complexes and subcomplexes. */

  /*! \ingroup plexGroup
    \brief Unit providing complexes and subcomplexes. */
  class plexUnit : 
    public mzr::unit
  {
  public:
    typedef cpx::plexQuery<mzrPlexSpecies,
			   mzrOmniPlex> plexQueryType;
  
  private:
    friend class mzrRecognizer;
    
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;

    // Maps %pairs of structural sites (i.e. particular sites on particular
    // mols) to sensitive points for dimerizations.
    //
    // This function is purely structural, unlike the mapping of %pairs of
    // site params to binding params.  That map is used to take the internal
    // state of molecules into account when computing binding parameters.
    // It's assumed now that binding parameters are determined by the shapes
    // (siteParams) assigned to the sites that are bound together.
    //
    // This is used e.g. in the constructor of plexFamily to connect the
    // plexFamily to the decomposition reactions that the plexFamily
    // participates in.
    //
    // Right now, it looks like every plexFamily has a reference to this map,
    // used whenever a new plexSpecies appears in the plexFamily.
    cpx::knownBindings<bnd::mzrMol, fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> > > bindingFeatures;

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
    std::set<mzrPlexFamily*> omniPlexFamilies;

    // Map from omniplexes (frequently parsed by the plexUnit for its client
    // units) to the nodes that were parsed to produce them (by the plexUnit
    // for its client units.  These are looked up later for the clients
    // using this map.
    std::map<xmlpp::Node*,
	     mzrOmniPlex*> nodeToOmni;

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
    utl::autoVector<fnd::query<cpx::omniStructureQueryArg<mzrPlex> > > structureQueries;

  public:
    plexUnit(mzr::moleculizer& rMoleculizer,
	     mzr::mzrUnit& refMzrUnit,
	     bnd::molUnit& refMolUnit);

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
    fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> >*
    addBindingFeature(bnd::mzrMol* pLeftMol,
		      int leftMolSiteSpec,
		      bnd::mzrMol* pRightMol,
		      int rightMolSiteSpec);

    /*! \brief Recover binding feature for %pair of structural sites.

    This is used to connect a plexFamily to its binding features.
    Doing so enables the plexFamily to notify the reaction families
    that are interested in bindings (only decompFam ?) when a new
    member species has been added to the plexFamily. */
    fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> >*
    findBindingFeature(bnd::mzrMol* pLeftMol,
		       int leftMolSiteSpec,
		       bnd::mzrMol* pRightMol,
		       int rightMolSiteSpec);


    // When an omniplex is parsed (possibly at the behest of another unit)
    // it should be added here.
    void
    addOmniPlex(mzrOmniPlex* pOmniPlex,
		xmlpp::Node* pParentNode)
      throw(utl::xcpt);

    // For use by plexFamily::connectToFeatures.  This allows one to determine
    // which omniplexes (i.e. distinguished subcomplexes) appear in a new
    // structural family of complex species.
    const std::set<mzrPlexFamily*>&
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
    mzrOmniPlex*
    getOmniForNode(xmlpp::Node* pParentNode) const;

    // Looks up the omniPlex parsed for a particular node (there is no
    // "omniPlex node," but an omniplex's nodes are all beneath the same
    // parent, which is the node argument here.)  Throws an "internal
    // exception" if no omniplex is found.
    mzrOmniPlex*
    mustGetOmniForNode(xmlpp::Node* pParentNode) const
      throw(utl::xcpt);

    mzrRecognizer recognize;

    void
    addPlexQuery(plexQueryType* pQuery)
    {
      plexQueries.push_back(pQuery);
    }

    void
    addStructureQuery(fnd::query<cpx::omniStructureQueryArg<mzrPlex> >* pQuery)
    {
      structureQueries.push_back(pQuery);
    }

//     bool
//     addPlexDumpable(const std::string& rDumpableName,
// 		    mzr::paramDumpable<plexSpecies>* pPlexDumpable)
//     {
//       return plexDumpables.addEntry(rDumpableName,
// 				    pPlexDumpable);
//     }

//     mzr::paramDumpable<plexSpecies>*
//     findPlexDumpable(const std::string& rDumpableName)
//     {
//       return plexDumpables.findEntry(rDumpableName);
//     }

//     bool
//     addOmniPlexDumpable(const std::string& rDumpableName,
// 			omniPlexDumpable* pDumpable)
//     {
//       return omniPlexDumpables.addEntry(rDumpableName,
// 					pDumpable);
//     }

//     omniPlexDumpable*
//     findOmniPlexDumpable(const std::string& rDumpableName)
//     {
//       return omniPlexDumpables.findEntry(rDumpableName);
//     }

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

#endif // PLEXUNIT_H
