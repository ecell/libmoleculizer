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
#include "mol/mol.hh"
#include "plex/plexFamily.hh"
#include "plex/plexFeature.hh"
#include "plex/plexUnit.hh"
#include "dimer/dimerUnit.hh"

namespace plx
{
  // Connects the plex family to its free site features; i.e. fills in
  // the freeSiteFeatures map.
  class plexFamily::doBehaviorizeSite
    : public std::unary_function<plexSiteSpec, void>
  {
    plexFamily& rFamily;
  public:
    doBehaviorizeSite(plexFamily& rPlexFamily) :
      rFamily(rPlexFamily)
    {}

    void
    operator()(const plexSiteSpec& rSpec) const
    {
      // Locate the default site parameter and the sensitiveFeature
      // for the "structural site".
      bnd::mol* pMol = rFamily.getParadigm().mols[rSpec.molNdx()];
      mzr::feature<plexSpecies, plexSiteSpec>& rSiteFeature
	= pMol->getSiteFeature(rSpec.siteNdx());

      // Install the site feature.
      rFamily.freeSiteFeatures[rSpec] = &rSiteFeature;
    }
  };

  // Install the binding feature for the given binding into the
  // binding feature map.
  void
  plexFamily::behaviorizeBinding(const plexBindingSpec& rSpec)
  {
    const plexBinding& rbinding = paradigm.bindings[rSpec];

    const plexSiteSpec& rLeftSpec = rbinding.leftSite();
    bnd::mol* pLeftMol = paradigm.mols[rLeftSpec.molNdx()];

    const plexSiteSpec& rRightSpec = rbinding.rightSite();
    bnd::mol* pRightMol = paradigm.mols[rRightSpec.molNdx()];

    // Try to look up the binding feature in the table.
    // findBindingFeature tries the "edge" in both directions.
    bindingFeature* pFeature
      = rPlexUnit.findBindingFeature(pLeftMol,
				     rLeftSpec.siteNdx(),
				     pRightMol,
				     rRightSpec.siteNdx());
    if(! pFeature)
      {
	// This means that the user specified a complex containing
	// a binding without there being a corresponding dimerization-gen.
	const bnd::bindingSite& rLeftSite
	  = pLeftMol->getSiteByNdx(rLeftSpec.siteNdx());
	const bnd::bindingSite& rRightSite
	  = pRightMol->getSiteByNdx(rRightSpec.siteNdx());

	throw noKineticConstsXcpt(pLeftMol->getName(),
				  rLeftSite.getName(),
				  pRightMol->getName(),
				  rRightSite.getName());
      }
    else
      {
	// Install the binding feature.
	bindingFeatures[rSpec] = pFeature;
      }
  }

  void
  plexFamily::behaviorizeMol(const plexMolSpec& rSpec)
  {
    bnd::mol* pMol = paradigm.mols[rSpec];

    // Install the sensitive feature.
    molFeatures[rSpec] = pMol;
  }

  class plexFamily::connectFamilyToSatisfiedOmni :
    public std::unary_function<omniPlex*, void>
  {
    plexFamily& rFamily;
    const plexIsoPair& rInjection;
  public:
    connectFamilyToSatisfiedOmni(plexFamily& rPlexFamily,
				 const plexIsoPair& rIsoPair) :
      rFamily(rPlexFamily),
      rInjection(rIsoPair)
    {}

    void
    operator()(omniPlex* pOmniPlex) const
    {
      const andOmniStructureQueries& rQuery = pOmniPlex->getStructureQuery();
      if(rQuery(rFamily.getParadigm(),
		rInjection))
	{
	  // Add the omniPlex's feature to the featureMap of rFamily.
	  rFamily.omniFeatures.addFeature(subPlexSpec(pOmniPlex,
						      rInjection),
					  pOmniPlex->getSubPlexFeature());
	}
    }
  };

  // This function both searches for subPlexes in the plexFamily's paradigm
  // and installs the ones that it finds.
  class plexFamily::doBehaviorizeOmni
    : public std::unary_function<plexFamily*, void>
  {
    plexFamily& rFamily;
  public:
    doBehaviorizeOmni(plexFamily& rPlexFamily) :
      rFamily(rPlexFamily)
    {}

    void
    operator()(plexFamily* pOmniFamily) const
    {
      // Is there an injection of the omni structure into this structure?
      plexIsoPair injection;
      if(reportIsoSearch(pOmniFamily->getParadigm(),
			 rFamily.getParadigm(),
			 injection).findInjection())
	{
	  // Attach this plex family to the features of those omniPlexes
	  // (associated to the omni family) whose structural queries
	  // are satisfied by the structure of rFamily.
	  std::for_each(pOmniFamily->omniPlexes.begin(),
			pOmniFamily->omniPlexes.end(),
			connectFamilyToSatisfiedOmni(rFamily,
						     injection));
	}
    }
  };

  // Connects this plexFamily to all its features, and thereby
  // to all the reactionFamilies that are sensitive to it by virtue
  // of connection to a structural feature.  For example, free
  // sites are features to which families of dimerization reactions
  // are sensitive.
  //
  // This should be done only after all the omniPlex families have
  // been put through passes 1 and 2.
  void
  plexFamily::connectToFeatures(void)
  {
    // Locate the free structural sites.
    std::vector<plexSiteSpec> freeSiteVector;
    paradigm.makeFreeSiteVector(freeSiteVector);

    // Install each free structural site's feature in the free-site
    // feature map.
    for_each(freeSiteVector.begin(),
	     freeSiteVector.end(),
	     doBehaviorizeSite(*this));

    // Install each binding's feature in the binding feature map.
    int bindingNdx = paradigm.bindings.size();
    while(0 < bindingNdx--) behaviorizeBinding((plexBindingSpec) bindingNdx);

    // Install each mol's feature in the mol feature map.
    int molNdx = paradigm.mols.size();
    while(0 < molNdx--) behaviorizeMol((plexMolSpec) molNdx);

    // Install the features for omniplexes.
    const std::set<plexFamily*>& rOmniFamilies
      = rPlexUnit.getOmniPlexFamilies();
    for_each(rOmniFamilies.begin(),
	     rOmniFamilies.end(),
	     doBehaviorizeOmni(*this));
  }
}
