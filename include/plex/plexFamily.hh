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

#ifndef PLEXFAMILY_H
#define PLEXFAMILY_H

/*! \defgroup plexStructGroup Structure
  \ingroup plexGroup
  \brief Structural equivalence, structural families of complexes. */

/*! \file plexFamily.hh
  \ingroup plexStructGroup
  \brief Defines plexFamily, a structural family of species of complexes. */

#include <vector>
#include <map>
#include <float.h>
#include <functional>
#include "mzr/util.hh"
#include "mzr/paramSpecies.hh"
#include "mzr/paramDumpable.hh"
#include "mzr/speciesFamily.hh"
#include "mzr/featureMap.hh"
#include "mol/molState.hh"
#include "plex/prm.hh"
#include "plex/plex.hh"
#include "plex/subPlexSpec.hh"
#include "plex/omniPlexFeature.hh"
#include "plex/alloSiteQuery.hh"
#include "plex/omniPlex.hh"
#include "plex/plexXcpt.hh"

namespace dimer
{
  class dimerUnit;
}

namespace plx
{
  class plexUnit;


  /*! \ingroup plexStructGroup

  \brief A structural family of species of complexes.

  The parameter that is used to classify complexes is the vector
  of molParams of the complex.  All other parameters of the
  complex are computed from these, together with the allostery
  properties of the structural family. */
  class plexFamily
    : public mzr::speciesFamily<std::vector<bnd::molParam>, plexSpecies>
  {
    //   There are two phases of manipulations to the plexFamilies that must
    //   be done before runtime:

    //   Pass 1:
    //     - construction
    //     - setting allosteric sites that are explicitly specified by
    //       "allo-site" commands.

    //   Pass 2:
    //     - connection to features("behaviorize"), including omniplex
    //     features.
    //     - setting allosteric sites from omniplexes, using the omniplex
    //     feature map to locate the relevant omniplexes.

    //   Passes 1 and 2 should be complete for omniplex families before
    //   pass 2 is done for other plex families.

  protected:
    // Some auxiliary function classes and adaptors.  These "friendly"
    // ones should become automatically friendly upon next compiler upgrade.
    //
    // The rest of them probably should not be listed here; just namespaced.
    class doBehaviorizeSite;
    friend class doBehaviorizeSite;

    class doBehaviorizeOmni;
    friend class doBehaviorizeOmni;

    class connectFamilyToSatisfiedOmni;
    friend class connectFamilyToSatisfiedOmni;
    
    class getDefaultMolParam;
    class injectAlloSiteIn; friend class injectAlloSiteIn;
    class applyOmniMods;
    class insertFreeSiteParam;
    class insertBindingParam;
    // Not sure, but this one looks obsolete.
    class dumpPop;
  
    // The first plex in this species that was seen.  This determines
    // the "official" ordering of the mols and bindings.
    plex paradigm;

    // The omniPlexes with structure given by the paradigm of this
    // plexFamily.  Putting these here makes it possible, after determining
    // that a new plexSpecies has structure with this plexFamily's structure
    // as a subcomplex, to apply the structural tests of just the omniPlexes
    // associated with this plexFamily's structure.
    //
    // These omniPlexes are memory managed by this plexFamily.
    mzr::autoVector<omniPlex> omniPlexes;

    // This feature should be watched by reactions that are
    // interested instances of this plex family that occur
    // as sub-plexes in instances of other plex families.
    //    omniPlexFeature subPlexFeature;

    // Map of site specs of free binding sites to corresponding features.
    // Sites themselves are now regarded like "symbols" that identify
    // a unique binding type: site pointers are used to look up binding
    // constants, so they amount to parameters.  The "allo-site" construct
    // modifies the default value of one of the parameters.
    mzr::featureMap<plexSpecies, plexSiteSpec> freeSiteFeatures;
  
    // Map of binding specs to corresponding features.
    mzr::featureMap<plexSpecies, plexBindingSpec> bindingFeatures;

    // Map of mol specs to corresponding features.
    mzr::featureMap<plexSpecies, plexMolSpec> molFeatures;

    // Map of subcomplex specs to corresponding features.
    mzr::featureMap<plexSpecies, subPlexSpec> omniFeatures;

    // In order to add state filtering to the allosteric-omni
    // and allosteric-plex constructs, we now give a "mapping"
    // from queries on plexParams to allosteric site maps.  Actually,
    // the mapping is a list of pairs, since we do need traversal
    // and don't need random access.
    //
    // When a new plex species is constructed, and we are figuring out its
    // binding site shapes in the "allostery" function, I now intend to run
    // through this list for each omniplex that has been recognized in the new
    // plex.  Whenever the state of the new plex satisfies a query, the
    // associateed "allosteric site map" from plexSiteSpecs to site shape
    // pointers is overlayed over the new species's siteParams.  When an
    // allosteric-omni or allosteric-plex specification is read, it may not
    // have any state specification in it, for compatibility with old models.
    // All plexFamilies must have the "passAll" map paired to what is
    // currently the alloSiteMap, now the default alloSiteMap.
    queryAllosteryList alloStateList;

    // Map giving the allosteric sites of this plexFamily.  These
    // sites may be bound or unbound.
    //    std::map<plexSiteSpec, bnd::siteParam> alloSiteMap;

    std::vector<mzr::paramDumpable<plexSpecies>*> familyDumpables;

    /////////////////// 
  
    // Attaches this plexFamily to the bindingFeatures of all its bindings.
    //
    // This enables the family to notify reaction families of new member
    // species of this plexFamily. For example, there is a decomposition
    // reaction family that handles decompositions of bindings between
    // these two mols at these two sites.  This allows each such
    // reaction family to construct new reactions for the new member
    // species.
    void
    behaviorizeBinding(const plexBindingSpec& rSpec);

    // Attaches this plexFamily to the molFeatures of all its mols.
    // See behaviorizeBinding above.
    void
    behaviorizeMol(const plexMolSpec& rSpec);

    // For modifying one site parameter in the database of default
    // parameters of all sites, both free and bound.
//     void
//     setAlloSite(const plexSiteSpec& rSpec,
// 		bnd::siteParam prm)
//     {
//       mzr::forceInsert(alloSiteMap,
// 		       std::pair<plexSiteSpec, bnd::siteParam>(rSpec, prm));
//     }

    // Why is this here?
    plexUnit& rPlexUnit;

  public:

    // Constructs a plexFamily for the combinatorial structure given
    // by rParadigm, which becomes the paradigm of the new plexFamily.
    //
    // This minimal constructor doesn't calculate the default parameter,
    // for example, and is used prior to looking for subplexes, etc.
    plexFamily(const plex& rParadigm,
	       plexUnit& refPlexUnit);

    // Notify the features, dumpables of a new species.
    void
    notify(plexSpecies* pSpecies)
    {
      std::cerr << "plexFamily "
		<< this
		<< " notifying for species "
		<< pSpecies
		<< "."
		<< std::endl;
      
      freeSiteFeatures.notifyNew(pSpecies);
      bindingFeatures.notifyNew(pSpecies);
      molFeatures.notifyNew(pSpecies);
      omniFeatures.notifyNew(pSpecies);

      std::for_each(familyDumpables.begin(),
		    familyDumpables.end(),
		    std::bind2nd(std::mem_fun(&mzr::paramDumpable<plexSpecies>::notify),
				 pSpecies));
    }

    // For accumulating all the plexSpecies in order to update them
    // by zero when regenerating reaction network.
    class accumulateOneSpecies : public
    std::unary_function<plexFamily::value_type, void>
    {
      std::vector<plexSpecies*>& rAllSpecies;
    public:
      accumulateOneSpecies(std::vector<plexSpecies*>& rAllPlexSpecies) :
	rAllSpecies(rAllPlexSpecies)
      {}

      void
      operator()(const argument_type& rEntry) const
      {
	rAllSpecies.push_back(rEntry.second);
      }
    };

    void
    accumulateSpecies(std::vector<plexSpecies*>& rAllSpecies)
    {
      std::for_each(begin(),
		    end(),
		    accumulateOneSpecies(rAllSpecies));
    }

    const plex&
    getParadigm(void) const
    {
      return paradigm;
    }

    plexSpecies*
    makeMember(const std::vector<bnd::molParam>& rMolParams);

    // This key routine fills in the gaps in plexParam, starting from
    // the molParams.  This routine controls how allosteric effects
    // depend on the structure of the complex.
    plexParam
    allostery(const std::vector<bnd::molParam>& rMolParams) const;

    // This returns a pointer, since it's destined to be inserted in
    // a featureMap...
//     omniPlexFeature*
//     getSubPlexFeature(void)
//     {
//       return &subPlexFeature;
//     }

    // Inserts an omniPlex, which should be in this plexFamily's structural
    // class, into this plexFamily's list of omniPlexes.  Later, a new
    // structural family that has this plexFamily's structure as a subcomplex
    // can run down this list looking for omniPlexes whose structural queries
    // are satisfied, and connect to their omniPlexFeatures.
    //
    // This plexFamily memory-manages its list of omniPlexes, too.
    //
    // This is used by the omniPlex constructor, which automatically
    // installs new omniPlexes into their plexFamilies.
    void
    addOmniPlex(omniPlex* pOmniPlex)
    {
      omniPlexes.push_back(pOmniPlex);
    }

    // This routine is probably obsoleted by the new treatment of omniPlexes.
    //
    // It's only used now in plexFamily code, I think, and as an accessor
    // of a private variable, it should go away.
    const queryAllosteryList&
    getAlloStateList(void) const
    {
      return alloStateList;
    }

    omniPlex*
    mustFindOmniForNode(xmlpp::Node* pParentNode) const
      throw(mzr::mzrXcpt);

    // This is so that one can set the allostery properties of this
    // plexFamily.
    void
    addAlloQueryAndMap(const andPlexQueries* pQuery,
		       const siteToShapeMap& rSiteToShapeMap)
    {
      alloStateList.addQueryAndMap(pQuery,
				   rSiteToShapeMap);
    }

    // Installs the allosteric sites given by the map.  This is used
    // to install allosteric sites given explicitly by the "allo-site"
    // command.
//     void
//     setAlloSites(const std::map<plexSiteSpec, bnd::siteParam>& rAlloSiteMap)
//     {
//       mzr::forceInsert(alloSiteMap,
// 		  rAlloSiteMap.begin(),
// 		  rAlloSiteMap.end());
//     }

    // Same as above, except remaps the keys using a plexMap.  This is
    // used to do the allosteric modifications dictated by omniPlexes
    // found in this complex, using the injection of the omniPlex into
    // this plex.
//     void
//     injectAlloSites(const std::map<plexSiteSpec, bnd::siteParam>& rAlloSiteMap,
// 		    const plexIsoPair& rInjection);

    // For getting the site-to-parameter map en masse.  This is used
    // to get the allosteric modifications dictated by an omniPlex
    // in order to install them in another plexFamily using the above
    // routine.
//     const std::map<plexSiteSpec, bnd::siteParam>&
//     getAlloSites(void)
//     {
//       return alloSiteMap;
//     }

    // Connects this plexFamily to all its features, and thereby
    // to all the reactionFamilies that are sensitive to it by virtue
    // of connection to a structural feature.  For example, free
    // sites are features to which families of dimerization reactions
    // are sensitive.
    //
    // This should be done only after all the omniPlexFamilies have
    // been put through passes 1 and 2.
    void
    connectToFeatures(void);

    // Makes the vector of default molParams.  This is used by plexUnit
    // to construct the default species of each named plex family, so
    // that the default species can be referred to with the name given
    // to the plex family.
    std::vector<bnd::molParam>
    makeDefaultMolParams(void) const;

    // Gives all the different ways in which the subcomplex with family
    // pOmniFamily occurs in complexes in this family.
    //
    // It could also be used as a predicate to simply determine
    // if a subcomplex relation holds between known species of complex.
    // For that, it's heavier weight than needed, but still probably
    // cheaper than doing a (redundant) isomorphism search.
    // Search the subcomplex features of this family for the given subcomplex.
    //
    // This added in support of the new version of nucleotide exchange reaction.
    //
    // It points to my stupidity in the featureMap template: it should be
    // templated on the actual type of the feature and make (or enforce) the
    // assumption that the actual feature type inherit from
    // feature<bearerSpecies, featureSpec>.  In this case, I get away with it,
    // because subPlexSpec actually includes the subcomplex's plexFamily*.
    std::vector<subPlexSpec>
    getOmniImbeddings(const plexFamily* pOmniFamily) const;

    // Add a query-based dumpable.
    void
    addDumpable(mzr::paramDumpable<plexSpecies>* pDumpable)
    {
      familyDumpables.push_back(pDumpable);
    }

    // Output routine.
    void
    insertSpecies(xmlpp::Element* pExplicitSpeciesElt,
		  double molarFactor) const
      throw(std::exception);
  };
}

#endif

