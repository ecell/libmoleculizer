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
#include "mol/siteShape.hh"
#include "mol/molState.hh"
#include "plex/plexFamily.hh"
#include "dimer/dimerUnit.hh"

namespace plx
{
  // Applies to the parameter for a new plexSpecies the allosteric
  // modifications that come from the omniplexes that have been recognized
  // in the species.
  class plexFamily::applyOmniMods :
    public std::unary_function
  <mzr::featureMap<plexSpecies, subPlexSpec>::value_type, void>
  {
    plexParam& rTargetParam;

  public:

    applyOmniMods(plexParam& rTargetPlexParam) :
      rTargetParam(rTargetPlexParam)
    {}

    void operator()
      (const mzr::featureMap<plexSpecies, subPlexSpec>::value_type& rEntry) const
    {
      const subPlexSpec& rSpec = rEntry.first;
      omniPlex* pOmni = rSpec.getOmni();

      // Does the species pass the omniPlex's state query?
      const andPlexQueries* pQuery = pOmni->getStateQuery();
      if(pQuery->applyTracked(rTargetParam,
			      rSpec))
	{
	  // Get the omniplex's site to shape map.
	  const siteToShapeMap& rOmniSiteToShapeMap
	    = pOmni->getSiteToShapeMap();

	  // Insert the allosteric site shapes.
	  rTargetParam.siteParams.setSiteShapes(rOmniSiteToShapeMap,
						rSpec);
	}
    }
  };

  /*! \ingroup plexStructGroup

  \brief Fills in rest of plexParam, given molParams, making allosteric
  modifications. */
  plexParam
  plexFamily::allostery(const std::vector<bnd::molParam>& rMolParams) const
  {
    plexParam result;

    // The molParams vector of the result is a copy of the given vector
    // of molParams.
    result.molParams = rMolParams;

    // Fill the database with the siteParams specified by the allostery
    // maps of the mols.
    for(int molNdx = 0;
	molNdx < (int) paradigm.mols.size();
	molNdx++)
      {
	// Look up the allosteric forms of the sites on the mol,
	// as determined from the mol's state.
	bnd::mol* pMol = paradigm.mols[molNdx];
	const bnd::molState* pMolState = rMolParams[molNdx];
	const std::vector<bnd::siteParam>& rSiteParams
	  = pMol->allostery(pMolState);

	// Install the allosteric siteParams into result.siteParams.
	for(int siteNdx = 0;
	    siteNdx < (int) pMol->getSiteCount();
	    siteNdx++)
	  {
	    result.siteParams.insert
	      (std::make_pair(plexSiteSpec(molNdx, siteNdx),
			      rSiteParams[siteNdx]));
	  }
      }

    // Modify the siteParams database as indicated by omniPlexes,
    // working "through" the injection.
    for_each(omniFeatures.begin(),
	     omniFeatures.end(),
	     applyOmniMods(result));

    // Modify siteParams database as indicated by this plexFamily.
    // Similar to what is done for each omniPlex above, but without the
    // permutation.
    getAlloStateList().setSatisfiedQuerySiteShapes(result);

    return result;
  }

  // This is for setting the allosteric sites that come from omniplexes.
  // First, it applies the omniplex's injection into this plexFamily
  // to the siteSpec to get a siteSpec valid in this plexFamily.  Then,
  // it installs the siteParam into this plexFamily's allostery map
  // "alloSiteMap" using the transformed siteSpec.
//   class plexFamily::injectAlloSiteIn :
//     public std::unary_function<std::map<plexSiteSpec, bnd::siteParam>::value_type, void>
//   {
//     plexFamily& rFamily;
//     const plexIsoPair& rInj;
//   public:
//     injectAlloSiteIn(plexFamily& rPlexFamily,
// 		     const plexIsoPair& rInjection) :
//       rFamily(rPlexFamily),
//       rInj(rInjection)
//     {}
//     void operator()
//       (const std::map<plexSiteSpec, bnd::siteParam>::value_type& rEntry) const
//     {
//       rFamily.setAlloSite(rInj.forward.applyToSiteSpec(rEntry.first),
// 			  rEntry.second);
//     }
//   };

  // Same as "setAlloSites", except remaps the keys using a plexMap.
  // This is used to do the allosteric modifications dictated by
  // omniPlexes found in this complex, using the injection of the
  // omniPlex into this plex.  This should be done in pass 2, after all
  // the omniplexes have been recognized.
//   void
//   plexFamily::injectAlloSites(const std::map<plexSiteSpec, bnd::siteParam>& rAlloSiteMap,
// 			      const plexIsoPair& rInjection)
//   {
//     for_each(rAlloSiteMap.begin(),
// 	     rAlloSiteMap.end(),
// 	     injectAlloSiteIn(*this,
// 			      rInjection));
//   }

  class plexFamily::getDefaultMolParam
    : public std::unary_function<bnd::mol*, bnd::molParam>
  {
  public:
    bnd::molParam
    operator()(bnd::mol* pMol) const
    {
      return pMol->getDefaultParam();
    }
  };

  std::vector<bnd::molParam>
  plexFamily::makeDefaultMolParams(void) const
  {
    std::vector<bnd::molParam> defaultMolParams;

    std::back_insert_iterator<std::vector<bnd::molParam> > outIter(defaultMolParams);

    transform(paradigm.mols.begin(),
	      paradigm.mols.end(),
	      outIter,
	      getDefaultMolParam());

    return defaultMolParams;
  }
}
