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

#ifndef GLOBALSPECIES_H
#define GLOBALSPECIES_H

#include "fnd/massive.hh"
#include "cpt/compartmentGraph.hh"
#include "cpt/compartmentSpecies.hh"
#include "fnd/notifier.hh"
#include "fnd/sensitive.hh"
#include "fnd/newSpeciesStimulus.hh"
#include "utl/string.hh"
#include "utl/gsl.hh"

namespace cpt
{
  // When the first molecule of an individual compartment species is created,
  // the individual compartment species notifies its globalSpecies, which
  // passes on the new-species notification.
  //
  // (For example, if this were a global species of complex, the global
  // species would pass along the notification for its compartment species,
  // along with the compartment index.  We expect all features of global
  // species to include the compartment index in their featurSpec.).
  //
  // This class is not copyable, since it manages its compartmentSpecies.
  class globalSpecies :
    public fnd::sensitive<fnd::newSpeciesStimulus<compartmentSpecies> >,
    public fnd::massive
  {
  protected:
    const compartmentGraph& rGraph;

    // Parallel to the compartments in the compartmentGraph.
    utl::autoVector<compartmentSpecies> compartmentSpeciesVector;

    // Parallel to the boundaries in the compartmentGraph.
    std::vector<double> diffusionRates;

    utl::gsl::autoGslMatrix diffusionMatrix;

    static int
    globalSpeciesCount;

  public:
    // This constructor might not ever be used.
    // 
    // The vector giving the compartment populations should be indexed
    // like the compartments in the compartment graph.
    //
    // The vector giving the diffusion rates should be indexed like the
    // boundaries in the compartment graph.
    globalSpecies(const compartmentGraph& rCompartmentGraph,
		  const std::vector<double>& rDiffusionRates,
		  const std::vector<int>& rCompartmentPops);

    // The vector giving the diffusion rates should be indexed like the
    // boundaries in the compartment graph.
    //
    // This could also be a static wrapper of the more flexible constructor
    // above.
    //
    // This constructor is used in parsing the species, where the default
    // population of 0 is taken, and another parsing pass initializes the
    // population.  This sounds like bad design; it comes from a way of doing
    // things imposed by plexUnit.  See "prepareToRun".
    globalSpecies(const compartmentGraph& rCompartmentGraph,
		  const std::vector<double>& rDiffusionRates,
		  int defaultCompartmentPop = 0);
    
    virtual
    ~globalSpecies(void)
    {}

    static int
    getGlobalSpeciesCount(void)
    {
      return globalSpeciesCount;
    }

    std::string
    getTag(void) const
    {
      return utl::stringify<const globalSpecies*>(this);
    }

    virtual
    std::string
    getName(void) const
    {
      return getTag();
    }

    const compartmentGraph&
    getCompartmentGraph(void) const
    {
      return rGraph;
    }

    // Response to first update of a subordinate compartment species,
    // ultimately supposed to lead to species/reaction generation.  This might
    // consist, for example of notifying, the plexFamily of a globalSpecies of
    // plexes.
    virtual
    void
    respond(const fnd::newSpeciesStimulus<compartmentSpecies>& rStim)
    {}

    void
    getConcVector(utl::gsl::autoGslVector& rConcVector) const;

    utl::gsl::autoGslVector
    getConcVector(void) const
    {
      utl::gsl::autoGslVector concVector(rGraph.compartments.size());
      getConcVector(concVector);
      return concVector;
    }

    // Gets the sum of the populations of the compartment species.
    int
    getTotalPop(void) const;

    // Gets the compartment species for a given compartment index.
    compartmentSpecies*
    getCompartmentSpecies(int compartmentNdx) const
    {
      return compartmentSpeciesVector[compartmentNdx];
    }

    // The target vector here must be of the correct size, which is, for
    // example, this->getCompartmentGraph().compartments.size().
    void
    getPopVector(std::vector<int>& rTargetVector) const
    {
      std::transform(compartmentSpeciesVector.begin(),
		     compartmentSpeciesVector.end(),
		     rTargetVector.begin(),
		     std::mem_fun(&compartmentSpecies::getPop));
    }

    std::vector<int>
    getPopVector(void) const
    {
      std::vector<int> result(compartmentSpeciesVector.size());
      getPopVector(result);
      return result;
    }

    void
    update(const std::vector<int>& compartmentDeltas,
	   fnd::sensitivityList<cptReaction>& rAffectedReactions,
	   int notifyDepth);

    void
    update(const std::map<int, int>& compartmentToDelta,
	   fnd::sensitivityList<cptReaction>& rAffectedReactions,
	   int notifyDepth);

    const utl::gsl::autoGslMatrix&
    getDiffusionMatrix(void) const
    {
      return diffusionMatrix;
    }

    void
    doDiffusion(utl::gsl::multinomialSampler& rSampler,
		double timeInterval,
		double epsilon,
		fnd::sensitivityList<cptReaction>& rAffectedReactions,
		int notifyDepth);

    virtual
    void
    respond(const compartmentSpecies* pUpdatedCptSpecies)
    {
      // Global stoch species won't do anything here, global
      // plex species will notifiy their plexFamily.
    }
  };
}

#endif // GLOBALSPECIES_H
