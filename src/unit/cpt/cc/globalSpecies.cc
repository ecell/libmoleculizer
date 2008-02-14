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

#include <functional>
#include "utl/forBoth.hh"
#include "cpt/globalSpecies.hh"

namespace cpt
{
  class insertBoundaryInMatrix :
    public std::unary_function<boundary, void>
  {
    utl::gsl::autoGslMatrix& rDiffusionMatrix;
  public:
    insertBoundaryInMatrix(utl::gsl::autoGslMatrix& refDiffusionMatrix) :
      rDiffusionMatrix(refDiffusionMatrix)
    {}

    void
    operator()(const boundary* pBoundary,
	       double diffusionRate) const
    {
      compartment* pFirstCpt = pBoundary->getFirstCpt();
      int firstNdx = pFirstCpt->getIndex();
      double firstVol = pFirstCpt->getVolume();

      compartment* pSecondCpt = pBoundary->getSecondCpt();
      int secondNdx = pSecondCpt->getIndex();
      double secondVol = pSecondCpt->getVolume();

      // This points out how we really only need one real parameter
      // to describe the boundary, the ratio of "area" to "thickness."
      // Comes in units of length.
      double boundaryRate
	= diffusionRate * pBoundary->getArea() / pBoundary->getThickness();

      // Should I be using the molar constant here, instead of the volume?
      // This depends on the units of the diffusion rate.
      rDiffusionMatrix.at(firstNdx,
			  firstNdx) -= boundaryRate / firstVol;
      rDiffusionMatrix.at(firstNdx,
			  secondNdx) += boundaryRate / secondVol;
      rDiffusionMatrix.at(secondNdx,
			  firstNdx) += boundaryRate / firstVol;
      rDiffusionMatrix.at(secondNdx,
			  secondNdx) -= boundaryRate / secondVol;
    }
  };

  globalSpecies::
  globalSpecies(const compartmentGraph& rCompartmentGraph,
		const std::vector<double>& rDiffusionRates,
		const std::vector<int>& rCompartmentPops) :
    rGraph(rCompartmentGraph),
    diffusionRates(rDiffusionRates),
    diffusionMatrix(rCompartmentGraph.compartments.size(),
		    rCompartmentGraph.compartments.size(),
		    true)
  {
    // Construct the compartment species.
    int compartmentCount = rCompartmentGraph.compartments.size();
    for(int compartmentNdx = 0;
	compartmentNdx < compartmentCount;
	++compartmentNdx)
      {
	compartmentSpeciesVector.push_back
	  (new compartmentSpecies(*this,
				  compartmentNdx,
				  rCompartmentPops[compartmentNdx]));
      }

    // Calculate the diffusion matrix, which "true" above caused to be
    // initialized to zero.
    utl::for_both(rGraph.boundaries.begin(),
		  rGraph.boundaries.end(),
		  diffusionRates.begin(),
		  insertBoundaryInMatrix(diffusionMatrix));

    // Bump the global species count for speciesCountDumpable.
    ++globalSpeciesCount;
  }

  globalSpecies::
  globalSpecies(const compartmentGraph& rCompartmentGraph,
		const std::vector<double>& rDiffusionRates,
		int defaultCompartmentPop) :
    rGraph(rCompartmentGraph),
    diffusionRates(rDiffusionRates),
    diffusionMatrix(rCompartmentGraph.compartments.size(),
		    rCompartmentGraph.compartments.size(),
		    true)
  {
    // Construct the compartment species.
    int compartmentCount = rCompartmentGraph.compartments.size();
    for(int compartmentNdx = 0;
	compartmentNdx < compartmentCount;
	++compartmentNdx)
      {
	compartmentSpeciesVector.push_back
	  (new compartmentSpecies(*this,
				  compartmentNdx,
				  defaultCompartmentPop));
      }

    // Calculate the diffusion matrix, which "true" above caused to be
    // initialized to zero.
    utl::for_both(rGraph.boundaries.begin(),
		  rGraph.boundaries.end(),
		  diffusionRates.begin(),
		  insertBoundaryInMatrix(diffusionMatrix));

    // Bump the global species count for speciesCountDumpable.
    ++globalSpeciesCount;
  }

  class accumulatePop :
    public std::unary_function<const compartmentSpecies*, void>
  {
    int& rTotal;
  public:
    accumulatePop(int& rTotalPop) :
      rTotal(rTotalPop)
    {}

    void
    operator()(const compartmentSpecies* pCompartmentSpecies) const
    {
      rTotal += pCompartmentSpecies->getPop();
    }
  };

  int
  globalSpecies::
  getTotalPop(void) const
  {
    int totalPop = 0;
    std::for_each(compartmentSpeciesVector.begin(),
		  compartmentSpeciesVector.end(),
		  accumulatePop(totalPop));
    return totalPop;
  }

  void
  globalSpecies::
  getConcVector(utl::gsl::autoGslVector& rConcVector) const
  {
    // Since interaction with gsl seems to be essentially index-based, there's
    // not much use in fooling with iterators/STL here.
    int compartmentNdx = rGraph.compartments.size();
    while(0 < compartmentNdx--)
      {
	rConcVector.at(compartmentNdx)
	  = compartmentSpeciesVector[compartmentNdx]->getConc();
      }
  }

  // This is for updating all the compartment species by given deltas, rather
  // than for resetting them to given populations. See
  // updateCompartmentSpeciesPop below.
  class updateSomeCompartmentSpecies :
    public std::unary_function<std::map<int, int>::value_type, void>
  {
    globalSpecies& rSpecies;
    fnd::sensitivityList<cptReaction>& rAffected;
    int depth;

  public:
    updateSomeCompartmentSpecies
    (globalSpecies& rGlobalSpecies,
     fnd::sensitivityList<cptReaction>& rAffectedReactions,
     int notifyDepth) :
      rSpecies(rGlobalSpecies),
      rAffected(rAffectedReactions)
    {}

    void
    operator()(argument_type cptNdxDeltaPair) const
    {
      // Get the compartment species using the given compartment index.
      int cptNdx = cptNdxDeltaPair.first;
      compartmentSpecies* pCptSpecies
	= rSpecies.getCompartmentSpecies(cptNdx);

      // Update its population by the given delta.
      int delta = cptNdxDeltaPair.second;
      pCptSpecies->update(delta,
			  rAffected,
			  depth);
    }
  };

  void
  globalSpecies::
  update(const std::map<int, int>& rCompartmentToDelta,
	 fnd::sensitivityList<cptReaction>& rAffectedReactions,
	 int notifyDepth)
  {
    std::for_each(rCompartmentToDelta.begin(),
		  rCompartmentToDelta.end(),
		  updateSomeCompartmentSpecies(*this,
					       rAffectedReactions,
					       notifyDepth));
  }

  // More missing stuff in our STL.  The single-argument version of
  // std::mem_fun can't accept a reference for its argument.  That's probably
  // something to do with partial specialization of templates.
  class updateCompartmentSpecies :
    public std::binary_function<compartmentSpecies*, int, void>
  {
    fnd::sensitivityList<cptReaction>& rAffected;
    int depth;

  public:
    updateCompartmentSpecies
    (fnd::sensitivityList<cptReaction>& rAffectedReactions,
     int notifyDepth) :
      rAffected(rAffectedReactions),
      depth(notifyDepth)
    {}
      
    void
    operator()(compartmentSpecies* pCSpecies,
	       int delta)
    {
      pCSpecies->update(delta,
			rAffected,
			depth);
    }
  };

  void
  globalSpecies::
  update(const std::vector<int>& rCompartmentDeltas,
	 fnd::sensitivityList<cptReaction>& rAffectedReactions,
	 int notifyDepth)
  {
    utl::for_both(compartmentSpeciesVector.begin(),
		  compartmentSpeciesVector.end(),
		  rCompartmentDeltas.begin(),
		  updateCompartmentSpecies(rAffectedReactions,
					   notifyDepth));
  }

  class setCompartmentSpeciesPop :
    public std::binary_function<int, compartmentSpecies*, void>
  {
    fnd::sensitivityList<cptReaction>& rAffected;
    int depth;
    
  public:
    setCompartmentSpeciesPop
    (fnd::sensitivityList<cptReaction>& rAffectedReactions,
     int notifyDepth) :
      rAffected(rAffectedReactions),
      depth(notifyDepth)
    {}
    
    void
    operator()(int newPop,
	       compartmentSpecies* pSpecies) const
    {
      pSpecies->update(newPop - pSpecies->getPop(),
		       rAffected,
		       depth);
    }
  };

  void
  globalSpecies::
  doDiffusion(utl::gsl::multinomialSampler& rSampler,
	      double timeInterval,
	      double epsilon,
	      fnd::sensitivityList<cptReaction>& rAffectedReactions,
	      int notifyDepth)
  {
    int compartmentCount = rGraph.compartments.size();
    
    // Convert the populations of the compartment species into
    // concentrations.
    utl::gsl::autoGslVector preConcentrations(compartmentCount);
    getConcVector(preConcentrations);

    // Exponentiate the diffusion matrix times the time interval.
    utl::gsl::autoGslMatrix solution(compartmentCount,
				     compartmentCount);
    diffusionMatrix.exp(timeInterval,
			epsilon,
			solution);

    // Examine the fates of the molecules that start out in each of the
    // compartments.
    // 
    // This amounts to extracting the columns, one by one, for multinomial
    // sampling.
    std::vector<int> postPopulations(compartmentCount, 0);
    int cptNdx = compartmentCount;
    while(0 < cptNdx--)
      {
	utl::gsl::autoGslVector cptFateVector(compartmentCount);
	solution.getColumnVector(cptNdx,
				 cptFateVector);

	// Sample multinomial to get counts of molecules coming
	// from compartment at cptNdx that land in all the other
	// compartments.
	std::vector<int> fatePops(compartmentCount);
	rSampler.sample(cptFateVector,
			getCompartmentSpecies(cptNdx)->getPop(),
			fatePops);

	std::transform(postPopulations.begin(),
		       postPopulations.end(),
		       fatePops.begin(),
		       postPopulations.begin(),
		       std::plus<int>());
      }

    // Reset the compartments' populations.
    // 
    // Note that, for now, this doesn't send any sort of update notification
    // to reactions.
    utl::for_both(postPopulations.begin(),
		  postPopulations.end(),
		  compartmentSpeciesVector.begin(),
		  setCompartmentSpeciesPop(rAffectedReactions,
					   notifyDepth));
  }

  int globalSpecies::globalSpeciesCount = 0;
}

