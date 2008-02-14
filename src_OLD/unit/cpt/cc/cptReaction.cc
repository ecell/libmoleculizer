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
#include "cpt/cptReaction.hh"
#include "cpt/respondReaction.hh"
#include "cpt/globalReaction.hh"
#include "cpt/cptUnit.hh"
#include "cpt/cptApp.hh"

namespace cpt
{
  cptReaction::
  cptReaction(globalReaction& refGlobalReaction,
	      int compartmentNdx) :
    fnd::basicReaction<compartmentSpecies>(),
    rGlobalReaction(refGlobalReaction),
    compartmentIndex(compartmentNdx)
  {
    // Sensitize reaction to changes in non-reactant state variables.
    //
    // In moleculizer, there was only the volume.  For now, I'm only going
    // to sensitize to the volume of the given compartment, leaving other
    // possible "global" state variables for later.  For the time being, it
    // looks like inter-compartment reactions (not necessarily a good idea
    // anyway) are on hold.

    // I'm not very happy with this way of getting to the compartmentGraph.
    compartment* pCompartment
      = rGlobalReaction.getCompartmentGraph().compartments[compartmentIndex];

    pCompartment->sensitizeReactionToVolume(*this);
  }

  compartment*
  cptReaction::
  getCompartment(void) const
  {
    return getGlobalReaction()
      .getCompartmentGraph()
      .compartments[getCompartmentIndex()];
  }

  // Now adds this reaction to the sensitivity list
  // of each reactant as it is added.
  void
  cptReaction::
  addReactant(compartmentSpecies* pSpecies,
	      int multiplicity)
  {
    fnd::basicReaction<compartmentSpecies>::addReactant(pSpecies,
						   multiplicity);
    pSpecies->addSensitive(this);
  }

  class combinationsForReactant :
    public std::unary_function<cptReaction::multMap::value_type, void>
  {
    double& rCombinations;
  public:
    combinationsForReactant(double& rCombinationCount) :
      rCombinations(rCombinationCount)
    {}

    void
    operator()(const argument_type& rSpeciesMultPair)
    {
      int speciesPop = rSpeciesMultPair.first->getPop();
      int multiplicity = rSpeciesMultPair.second;

      // Note that if speciesPop drops to 0 during this operation,
      // then rCombinations becomes 0.
      while(0 < multiplicity--) rCombinations *= speciesPop--;
    }
  };

  // Returns the current propensity.
  double
  cptReaction::
  propensity(void) const
  {
    double combinations = 1.0;
    std::for_each(reactants.begin(),
		  reactants.end(),
		  combinationsForReactant(combinations));
    return
      combinations
      * getRate()
      / std::pow(getCompartment()->getVolume() * fnd::avogadrosNumber,
	    getArity() - 1);
  }

  class updateSpecies :
    public std::unary_function<cptReaction::multMap::value_type, void>
  {
    fnd::sensitivityList<cptReaction>& rSensitives;
    int depth;

  public:
    updateSpecies(fnd::sensitivityList<cptReaction>& rSensitiveReactions,
		  int generateDepth) :
      rSensitives(rSensitiveReactions),
      depth(generateDepth)
    {}

    void
    operator()(const argument_type& rSpeciesDeltaPair) const
    {
      compartmentSpecies* pSpecies = rSpeciesDeltaPair.first;
      int delta = rSpeciesDeltaPair.second;

      pSpecies->update(delta,
		       rSensitives,
		       depth);
    }
  };

  fnd::eventResult
  cptReaction::
  happen(cptApp& rApp)
    throw(std::exception)
  {
    // Adjust the populations of reactants and products using the deltas map.
    //
    // Accumulate the reactions that should respond to these changes in
    // population in a set, for the time being, so that no reaction is called
    // on to update its propensity more than once.
    fnd::sensitivityList<cptReaction> sensitiveReactions;
    std::for_each(deltas.begin(),
		  deltas.end(),
		  updateSpecies(sensitiveReactions,
				rApp.getCptUnit().getGenerateDepth()));

    // Respond all the sensitive reactions.
    std::for_each(sensitiveReactions.begin(),
		  sensitiveReactions.end(),
		  respondReaction(rApp.getPropensities()));

    // Move this reaction to the front of the distribution.  This is a cheap
    // way of keeping high-propensity reactions near the front of the
    // distribution; no kind of actual sorting is necessary.
    iDistroEntry = rApp.getPropensities().moveToFront(iDistroEntry);

    // Bump the reaction event count.
    ++eventCount;

    return fnd::go;
  }

  void
  cptReaction::
  respond(const reactionStim& rStim)
  {
    // Calculate the change in propensity.
    double delta = propensity() - iDistroEntry->first;

    // Update the propensity entry for this reaction.  Trying to get the sum
    // of all the propensities in the list always to match totalPropensity.
    iDistroEntry->first += delta;

    // Update the total propensity.
    rStim
      .rDistro
      .updateTotalPropensity(delta);
  }

  int cptReaction::eventCount = 0;
}
