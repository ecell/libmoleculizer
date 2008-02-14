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

#include "utl/gsl.hh"
#include "mzr/mzrSpecies.hh"
#include "mzr/mzrReaction.hh"
#include "mzr/respondReaction.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/dumpUtils.hh"
#include "mzr/unitsMgr.hh"
#include "mzr/moleculizer.hh"

namespace mzr
{
  void
  mzrReaction::
  addReactant(mzrSpecies* pSpecies,
	      int multiplicity)
  {
    fnd::gillspReaction<mzrSpecies>::addReactant(pSpecies,
						 multiplicity);
    pSpecies->addSensitive(this);
  }

  void
  mzrReaction::
  respond(const mzrReactionStimulus& rStimulus)
  {
    moleculizer& rMzr = rStimulus.rMzr;
    eventQueue& rQueue = rMzr.eventQ;
    mzrUnit& rMzrUnit = *(rMzr.pUserUnits->pMzrUnit);
    utl::gsl::autoGslRng& rRng = rMzrUnit.rng;
    
    // Now trying to avoid scheduling events at time "never", as it appears
    // that inserting events in the event queue is the single biggest task
    // in this simulator.

    // Now that most things use volume (since other things use compartment
    // volumes, instead of the unique volume, I should probably make the
    // volume directly available in mzrUnit.  An alternative would be for
    // compartments to maintain a molarFactor (i.e. avogadrosNumber * volume).
    double myPropensity = propensity(rMzrUnit.getMolarFactor().getVolume());

    if(0.0 < myPropensity)
      {
	if((lastPropensity <= 0.0)
	   || (myPropensity > lastPropensity * highSensitive)
	   || (myPropensity < lastPropensity * lowSensitive))
	  {
	    lastPropensity = myPropensity;

	    double nextReactionTime
	      = rQueue.getSimTime()
	      + utl::gsl::exponentialSampler(rRng).sample(1.0 / myPropensity);

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

  class updateOneSpecies :
    public std::unary_function<fnd::basicReaction<mzrSpecies>::multMap::value_type, void>
  {
    fnd::sensitivityList<mzrReaction>& rAffected;
    int depth;
    
  public:
    updateOneSpecies(fnd::sensitivityList<mzrReaction>& rAffectedReactions,
		     int generateDepth) :
      rAffected(rAffectedReactions),
      depth(generateDepth)
    {}

    void
    operator()(const argument_type& rSpeciesMultPair) const
    {
      mzrSpecies* pSpecies = rSpeciesMultPair.first;
      int multiplicity = rSpeciesMultPair.second;

      pSpecies->update(multiplicity,
		       rAffected,
		       depth);
    }
  };

  // What a reaction event does, given what the reaction does.
  fnd::eventResult
  mzrReaction::
  happen(moleculizer& rMolzer)
    throw(std::exception)
  {
    // Reactions whose reactants are affected by this reaction
    // are accumulated in this set by doReaction.
    fnd::sensitivityList<mzrReaction> affectedReactions;

    // This is a leftover of the hack for "no reactant" arrows.  Since
    // they're not sensitive to anything, they have to reschedule
    // themselves.
    respondReaction doRespondReaction(rMolzer);
    if(reactants.size() == 0)
      {
	// This is a hacky special case, and that's why I don't like it.
	lastPropensity = -1.0;
	doRespondReaction(this);
      }

    // Do the reaction deltas.
    int generateDepth
      = rMolzer.pUserUnits->pMzrUnit->getGenerateDepth();
    std::for_each(deltas.begin(),
		  deltas.end(),
		  updateOneSpecies(affectedReactions,
				   generateDepth));

    // The reaction event was descheduled before being executed.
    // This should probably be a part of descheduling.
    //
    // This is bad bad bad.  It was connected with the "tolerance"
    // acceleration, which has left such spoor all over the place.
    // 
    // Are there any other times that a reaction is descheduled,
    // besides this and when the reaction is executed?
    lastPropensity = -1.0;

    // Reschedule the reactions whose reactants were affected by this
    // reaction. 
    affectedReactions.forEachSensitive(doRespondReaction);

    // For the dumpable that gives the total number of reactions so far.
    reactionEventCount++;

    return fnd::go;
  }

  // Following are connected with generating state dump output.

  class insertReactantSpeciesRef :
    public std::unary_function<std::map<mzrSpecies*, int>::value_type, void>
  {
    xmlpp::Element* pReactionElt;
  public:
    insertReactantSpeciesRef(xmlpp::Element* pReactionElement) :
      pReactionElt(pReactionElement)
    {
    }

    void
    operator()(const argument_type& rSpeciesMultPair) const
      throw(std::exception)
    {
      // Insert element for reaction reactant species.
      xmlpp::Element* pReactantSpeciesRefElt
	= pReactionElt->add_child(eltName::taggedSubstrate);

      // Add the name or tag of the substrate species as attribute.
      pReactantSpeciesRefElt
	->set_attribute(eltName::taggedSubstrate_tagAttr,
			rSpeciesMultPair.first->getTag());

      // Add multiplicity of substrate as attribute.
      pReactantSpeciesRefElt
	->set_attribute(eltName::taggedSubstrate_multAttr,
			utl::stringify<int>(rSpeciesMultPair.second));
    }
  };

  // It looks like I'm going back to emitting reaction substrates and
  // products, rather than substrates and deltas.  This is mainly for the
  // "immovable object," SBML.
  class insertProductSpeciesRef : public
  std::unary_function<std::map<mzrSpecies*, int>::value_type, void>
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
			utl::stringify<int>(rSpeciesMultPair.second));
    }
  };

  xmlpp::Element*
  mzrReaction::insertElt(xmlpp::Element* pParentElt) const 
    throw(std::exception)
  {
    xmlpp::Element* pReactionElt
      = pParentElt->add_child(eltName::tagReaction);
      
    std::for_each(reactants.begin(),
		  reactants.end(),
		  insertReactantSpeciesRef(pReactionElt));

    std::for_each(products.begin(),
		  products.end(),
		  insertProductSpeciesRef(pReactionElt));

    // Additional scientific notation here for use by ECell.
    mzr::addDoubleParamChild(pReactionElt,
			     eltName::rate,
			     eltName::rate_valueAttr,
			     rate);
    return pReactionElt;
  }

  // Reaction related dumpable globals.
  int mzrReaction::reactionCount = 0;
  int mzrReaction::reactionEventCount = 0;

  // Globals connected with the tolerance optimization.
  double mzrReaction::lowSensitive = 1.0;
  double mzrReaction::highSensitive = 1.0;
}
