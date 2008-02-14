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

#include <cmath>
#include <limits>
#include "utl/forBoth.hh"
#include "clx/cptPlexFamily.hh"
#include "clx/cptOmniPlex.hh"

namespace clx
{
  cptPlexFamily::
  cptPlexFamily(const cptPlex& rParadigm,
		cpx::knownBindings<cml::cptMol, fnd::feature<cpx::cxBinding<cptPlexSpecies, cptPlexFamily> > >& refKnownBindings,
		std::set<cptPlexFamily*>& refOmniplexFamilies,
		cpt::cptUnit& refCptUnit,
		const cpt::compartmentGraph& rCompartmentGraph) :
    cpx::plexFamily<cml::cptMol,
		    cptPlex,
		    cptPlexSpecies,
		    cptPlexFamily,
		    cptOmniPlex>(rParadigm,
				 refKnownBindings,
				 refOmniplexFamilies),
    rGraph(rCompartmentGraph),
    rCptUnit(refCptUnit)
  {}

  class checkMolBoundaryRate :
    public std::binary_function<const cml::cptMol*, cpx::molParam, void>
  {
    // To be the minimum diffusion rate for this boundary,
    // minimum being taken over mols.
    double& rMinRate;

    // The weight of the mol having the minimum diffusion rate over this
    // boundary.
    double& rMinRateWeight;

    // The total mass of the mols in their given states.  This is the
    // same formula as plexSpecies use to compute their mass, but obviously
    // not using the same code!
    double& rTotWeight;

    // The index of this boundary.
    int boundaryNdx;

  public:
    checkMolBoundaryRate(double& refMinRate,
			 double& refMinRateWeight,
			 double& refTotWeight,
			 int boundaryIndex) :
      rMinRate(refMinRate),
      rMinRateWeight(refMinRateWeight),
      rTotWeight(refTotWeight),
      boundaryNdx(boundaryIndex)
    {}

    void
    operator()(const cml::cptMol* pCptMol,
	       cpx::molParam param) const
    {
      double molWeight = param->getMolWeight();

      rTotWeight += molWeight;
      
      double molBoundaryRate = pCptMol->getBoundaryRate(boundaryNdx);
      if(rMinRate > molBoundaryRate)
	{
	  rMinRate = molBoundaryRate;
	  rMinRateWeight = molWeight;
	}
    }
  };

  cptPlexSpecies*
  cptPlexFamily::
  constructSpecies(const cpx::siteToShapeMap& rSiteParams,
		   const std::vector<cpx::molParam>& rMolParams)
  {

    // Clearly, I screwed up compartmentGraph with respect to this indexing
    // of boundaries and compartments stuff.
    int boundaryCount = rGraph.boundaries.size();
    std::vector<double> boundaryRates(boundaryCount);
    for(int boundaryNdx = 0;
	boundaryNdx < boundaryCount;
	++boundaryNdx)
      {
	// This is clearly sub-optimal.  It would be much better if diffusion
	// rates, like site shapes, depended on the mol's modification state.
	//
	// To accomplish this, alloMol would have to be cusomized to have
	// diffusion rates in the values of its "alloMap," which now
	// just gives the binding site shapes for allosteric states.
	double minRate = std::numeric_limits<double>::max();
	double massOfMin = -1.0;
	double totalMass = 0.0;
	utl::for_both(paradigm.mols.begin(),
		      paradigm.mols.end(),
		      rMolParams.begin(),
		      checkMolBoundaryRate(minRate,
					   massOfMin,
					   totalMass,
					   boundaryNdx));

	double adjustedRate =  minRate * std::sqrt(massOfMin / totalMass);

	boundaryRates[boundaryNdx] = adjustedRate;
      }

    cptPlexSpecies* pNewSpecies 
      = new cptPlexSpecies(*this,
			   rSiteParams,
			   rMolParams,
			   rGraph,
			   boundaryRates);

    rCptUnit.addSpecies(pNewSpecies);
    
    return pNewSpecies;
  }

  class insertCptPlexSpecies :
    public std::unary_function<cptPlexFamily::value_type, void>
  {
    xmlpp::Element* pExplicitSpeciesElt;
    
  public:
    insertCptPlexSpecies(xmlpp::Element* pExplicitSpeciesElement) :
      pExplicitSpeciesElt(pExplicitSpeciesElement)
    {}

    void
    operator()(const argument_type& rPlexFamilyEntry) const
      throw(std::exception)
    {
      cptPlexSpecies* pSpecies = rPlexFamilyEntry.second;
      pSpecies->insertElt(pExplicitSpeciesElt);
    }
  };

  void
  cptPlexFamily::
  insertSpecies(xmlpp::Element* pExplicitSpeciesElt) const
    throw(std::exception)
  {
    std::for_each(begin(),
		  end(),
		  insertCptPlexSpecies(pExplicitSpeciesElt));
  }
}
