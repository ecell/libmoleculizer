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

#include <iostream>
#include <algorithm>
#include "mzr/species.hh"
#include "mzr/reaction.hh"
#include "sampleDist/sampleDist.hh"
#include "mzr/moleculizer.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/dumpUtils.hh"

namespace mzr
{
  void
  reaction::schedule(eventQueue& rQueue,
		     mzrUnit& rMzrUnit)
  {
    // Now trying to avoid scheduling events at time "never", as it appears
    // that inserting events in the event queue is the single biggest task
    // in this simulator.

    double myPropensity = propensity(rMzrUnit);

    if(0.0 < myPropensity)
      {
	if((lastPropensity <= 0.0)
	   || (myPropensity > lastPropensity * highSensitive)
	   || (myPropensity < lastPropensity * lowSensitive))
	  {
	    // Now we only count it as an activation if the event
	    // really reschedules.
	    lastPropensity = myPropensity;
	    reaction::activationCount++;

	    double delta = sampleDist::sampleExp() / myPropensity;
	    double nextReactionTime = rQueue.getSimTime() + delta;

	    rQueue.scheduleEvent(this, nextReactionTime);
	  }
      }
    else
      {
	// It would probably be better to make the reinitialization
	// of lastPropensity to -1.0 to be a part of descheduling.
	// For a first try, putting it here and where the reaction executes.
	//
	// Are there any other times that a reaction is descheduled,
	// besides this and when the reaction is executed?
	lastPropensity = -1.0;
	rQueue.descheduleEvent(this);
      }
  }

  // This function allows the parser to add substrates to the reaction
  // subsequent to its creation, avoiding the problem mentioned above
  // for the propensity in the constructor.
  void
  reaction::addSubstrate(species* pSpecies,
			 int multiplicity)
  {
    // Check to see if the species already has an entry in the substrates.
    // 
    // Note this can be done more efficiently with correct use of
    // insert: attempt the insert with just the bare multiplicity.  If
    // that fails, use the interator returned by insert to adjust the
    // multiplicity that is already in place.  This avoids searching
    // once for find, then a second time for insert.
    std::map<species*, int>::iterator iSpeciesMultPair
      = substrates.find(pSpecies);
    if(substrates.end() == iSpeciesMultPair)
      {
	// This is a new species, so we insert it with the given
	// (unchecked) multiplicity.
	substrates.insert(std::make_pair(pSpecies,
					 multiplicity));

	// Since this is a new substrate species, we add this reaction
	// to the substrate species's vector of sensitive reactions.
	pSpecies->addSensitiveReaction(this);
      }
    else
      {
	// This species has already occurred as a substrate.
	// This would be the case if somebody wrote an equation
	// like A + A -> AA @ 550.0.
	iSpeciesMultPair->second += multiplicity;
      }

    // Check to see if the species already has an entry in the
    // reactionDeltas; then create or update the entry.
    // 
    // Note this can be made more efficient with correct use of insert,
    // as mentioned above.
    std::map<species*, int>::iterator iSpeciesDeltaPair
      = reactionDeltas.find(pSpecies);
    if(reactionDeltas.end() == iSpeciesDeltaPair)
      {
	reactionDeltas.insert(std::make_pair(pSpecies,
					     -multiplicity));
      }
    else
      {
	iSpeciesDeltaPair->second -= multiplicity;
      }

    // Increase the arity by the multiplicity.
    theArity += multiplicity;
  }

  class countSubstrate :
    public std::unary_function<std::pair<species*, int>, void>
  {
    int& rCount;
  public:
    countSubstrate(int& rSubstrateCount) :
      rCount(rSubstrateCount)
    {}

    void
    operator()(const std::pair<species*, int>& rEntry) const
    {
      rCount += rEntry.second;
    }
  };

  int
  reaction::getSubstrateCount(void)
  {
    int totCount = 0;
    for_each(substrates.begin(),
	     substrates.end(),
	     countSubstrate(totCount));
    return totCount;
  }

  // This function allows the parser to add products to the reaction
  // subsequent to its creation, avoid the problem mentioned above
  // for the propensity in the constructor.
  void
  reaction::addProduct(species* pSpecies,
		       int multiplicity)
  {
    // Create or update an entry for the product species in products,
    // which is mainly used for dumping.
    std::pair<std::map<species*, int>::iterator, bool>
      insertProductResult = products.insert(std::make_pair(pSpecies,
							   multiplicity));
    // If the insert failed, there is already an entry in the products
    // map, so we increment the entry multiplicity.
    if(! insertProductResult.second)
      {
	insertProductResult.first->second += multiplicity;
      }

    // Create or update an entry for the product species in reactionDeltas.
    std::pair<std::map<species*, int>::iterator, bool>
      insertDeltaResult = reactionDeltas.insert(std::make_pair(pSpecies,
							       multiplicity));
    if(! insertDeltaResult.second)
      {
	insertDeltaResult.first->second += multiplicity;
      }
  }

  /*! \ingroup stochGroup

  \brief Auxiliary function class for reaction propensity calculation.

  This multiplies the number of reactive combinations of molecules
  by factor coming from one substrate's population p and multiplicity m
  in the reaction: p(p - 1)(p - 2)(p - m + 1).

  Note that this is NOT the actual number of combinations; rather it
  is the number of combinations times m!.*/
  class combinationsForSubstrate :
    public std::unary_function<std::pair<species*, int>, void>
  {
    double& rCombinations;
  public:
    combinationsForSubstrate(double& rCombinationCount) :
      rCombinations(rCombinationCount)
    {}
    void operator()(const std::pair<species*, int>& rSpeciesMult)
    {
      int speciesPop = rSpeciesMult.first->getPop();

      int multiplicity = rSpeciesMult.second;
      // Note that if speciesPop drops to 0 during this operation,
      // then rCombinations becomes 0.
      while(0 < multiplicity--) rCombinations *= speciesPop--;
    }
  };

  double
  reaction::propensity(mzrUnit& rMzrUnit) const
  {
    // Figure out the number of combinations of substrate molecules.
    // 
    // Note that this will get the correct number of combinations (1)
    // for a "no-substrate" reaction.
    double combinations = 1.0;
    for_each(substrates.begin(),
	     substrates.end(),
	     combinationsForSubstrate(combinations));

    // Here the factors of m! that were multiplied into the
    // number of combinations compensate for using the
    // ordinary determinisitc reaction rate.
    return combinations * rate
      / pow(rMzrUnit.getMolarFactor().getFactor(),
	    theArity - 1);
  }

  class updateOneSpecies :
    public std::unary_function<std::pair<species*, int>, void>
  {
    std::set<reaction*>& rAffected;
  public:
    updateOneSpecies(std::set<reaction*>& rAffectedReactions) :
      rAffected(rAffectedReactions)
    {}

    void
    operator()(const std::pair<species*, int>& rSpeciesDelta) const
    {
      rSpeciesDelta.first->update(rAffected,
				  rSpeciesDelta.second);
    }
  };

  // What a reaction event does, given what the reaction does.
  event::eventResult
  reaction::doEvent(moleculizer& rMolzer,
		    mzrUnit& rMzrUnit)
  {
    // Reactions whose substrates are affected by this reaction
    // are accumulated in this set by doReaction.
    std::set<reaction*> affectedReactions;

    // This is a leftover of the hack for "no substrate" arrows.  Since
    // they're not sensitive to anything, they have ro reschedule
    // themselves.  This may need to change now that sort of reaction
    // has taken over the universe.
    if(substrates.size() == 0)
      {
	// This is a hacky special case, and that's why I don't like it.
	lastPropensity = -1.0;
	schedule(rMolzer.eventQ,
		 rMzrUnit);
      }

    // Do the reaction deltas.
    for_each(reactionDeltas.begin(),
	     reactionDeltas.end(),
	     updateOneSpecies(affectedReactions));

    // The reaction event was descheduled before being executed.
    // This should probably be a part of descheduling.
    //
    // This is bad bad bad.  It was connected with the "tolerance"
    // acceleration, which has left such spoor all over the place.
    // 
    // Are there any other times that a reaction is descheduled,
    // besides this and when the reaction is executed?
    lastPropensity = -1.0;

    // Reschedule the reactions whose substrates were affected by this
    // reaction. 
    for_each(affectedReactions.begin(),
	     affectedReactions.end(),
	     scheduleReaction(rMolzer.eventQ,
			      rMzrUnit));

    // For the dumpable that gives the total number of reactions so far.
    reactionEventCount++;

    return event::go;
  }

  class sensitizeToOneSubstrate :
    public std::unary_function<std::pair<species* const, int>, void>
  {
    reaction* pReaction;
  public:
    sensitizeToOneSubstrate(reaction* ptrReaction) :
      pReaction(ptrReaction)
    {}
    void operator()(const std::pair<species* const, int>& speciesMultPair) const
    {
      species* const pSpecies = speciesMultPair.first;
      pSpecies->addSensitiveReaction(pReaction);
    }
  };

  // This base method should be invoked by all descendant classes
  // of reaction to sensitize themselves to globals, like volume
  // and temperature.
  void
  reaction::sensitizeToSubstrates(mzrUnit& rMzrUnit)
  {
    // Sensitize to globals, like volume.  I note that all reactions are
    // sensitized to this, but unary reactions don't really need to be.
    rMzrUnit.getMolarFactor().addSensitiveReaction(this);

    // Sensitize to reaction substrates.
    for_each(substrates.begin(),
	     substrates.end(),
	     sensitizeToOneSubstrate(this));
  }

  class insertSubstrateSpeciesRef :
    public std::unary_function<std::map<species*, int>::value_type, void>
  {
    xmlpp::Element* pReactionElt;
  public:
    insertSubstrateSpeciesRef(xmlpp::Element* pReactionElement) :
      pReactionElt(pReactionElement)
    {
    }

    void
    operator()(const argument_type& rSpeciesMultPair) const
      throw(std::exception)
    {
      // Insert element for reaction substrate species.
      xmlpp::Element* pSubstrateSpeciesRefElt
	= pReactionElt->add_child(eltName::taggedSubstrate);

      // Add the name or tag of the substrate species as attribute.
      pSubstrateSpeciesRefElt
	->set_attribute(eltName::taggedSubstrate_tagAttr,
			rSpeciesMultPair.first->getTag());

      // Add multiplicity of substrate as attribute.
      pSubstrateSpeciesRefElt
	->set_attribute(eltName::taggedSubstrate_multAttr,
			domUtils::stringify<int>(rSpeciesMultPair.second));
    }
  };

//   class insertDeltaSpeciesRef :
//     public std::unary_function<std::map<species*, int>::value_type, void>
//   {
//     xmlpp::Element* pReactionElt;
//   public:
//     insertDeltaSpeciesRef(xmlpp::Element* pReactionElement) :
//       pReactionElt(pReactionElement)
//     {}

//     void
//     operator()(const argument_type& rSpeciesMultPair) const
//       throw(std::exception)
//     {
//       // Insert element for reaction delta
//       xmlpp::Element* pDeltaSpeciesRefElt
// 	= pReactionElt->add_child(eltName::taggedDelta);

//       // Add species tag as attribute.
//       pDeltaSpeciesRefElt
// 	->set_attribute(eltName::taggedDelta_tagAttr,
// 			rSpeciesMultPair.first->getTag());

//       // Add delta as attribute.
//       pDeltaSpeciesRefElt
// 	->set_attribute(eltName::taggedDelta_deltaAttr,
// 			domUtils::stringify<int>(rSpeciesMultPair.second));
//     }
//   };

  // It looks like I'm going back to emitting reaction substrates and
  // products, rather than substrates and deltas.  This is mainly for the
  // "immovable object," SBML.
  class insertProductSpeciesRef : public
  std::unary_function<std::map<species*, int>::value_type, void>
  {
    xmlpp::Element* pReactionElt;
  public:
    insertProductSpeciesRef(xmlpp::Element* pReactionElement) :
      pReactionElt(pReactionElement)
    {}

    void
    operator()(const argument_type& rSpeciesMultPair) const
      throw(std::exception)
    {
      // Insert element for reaction product.
      xmlpp::Element* pProductSpeciesRefElt
	= pReactionElt->add_child(eltName::taggedProduct);

      // Add species tag as attribute.
      pProductSpeciesRefElt
	->set_attribute(eltName::taggedProduct_tagAttr,
			rSpeciesMultPair.first->getTag());

      // Add multiplicity as attribute.
      pProductSpeciesRefElt
	->set_attribute(eltName::taggedProduct_multAttr,
			domUtils::stringify<int>(rSpeciesMultPair.second));
    }
  };

  void
  reaction::insertElt(xmlpp::Element* pParentElt) const throw(std::exception)
  {
    xmlpp::Element* pReactionElt
      = pParentElt->add_child(eltName::tagReaction);
      
    std::for_each(substrates.begin(),
		  substrates.end(),
		  insertSubstrateSpeciesRef(pReactionElt));

    std::for_each(products.begin(),
		  products.end(),
		  insertProductSpeciesRef(pReactionElt));

    mzr::addDoubleParamChild(pReactionElt,
			     eltName::rate,
			     eltName::rate_valueAttr,
			     rate);
    
//     xmlpp::Element* pRateElt
//       = pReactionElt->add_child(eltName::rate);

//     pRateElt->set_attribute(eltName::rate_valueAttr,
// 			    domUtils::stringify<double>(rate));
  }

  int reaction::reactionCount = 0;
  double reaction::activationCount = 0.0;
  int reaction::reactionEventCount = 0;
  double reaction::lowSensitive = 1.0;
  double reaction::highSensitive = 1.0;
}
