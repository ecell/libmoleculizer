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

#ifndef REACTION_H
#define REACTION_H

/*! \defgroup reactionGroup Reactions
  \ingroup mzrGroup
  \brief Reactions and their families. */

/*! \file reaction.hh
  \ingroup reactionGroup
  \brief Defines basic reaction classes and scheduling functions. */

#include <set>
#include <list>
#include <algorithm>
#include "mzr/event.hh"

namespace mzr
{
  /*! \ingroup reactionGroup
    \brief Base class of all reactions. */
  class reaction :
    public event
  {
    static double lowSensitive;
    static double highSensitive;

    // For printing purposes.
    friend class reactionElement;

    // The substrate species, with multiplicities.  We need these to
    // calculate propensities.
    std::map<species*, int> substrates;

    // The product species with multiplicities.  Technically, we could
    // reconstruct these from the reactionDeltas and the substrates, but
    // keeping them in a separate pit makes it easier to dump by tags
    // (i.e. dump a network of Stochastirator reactions.)
    std::map<species*, int> products;

    // Rather than treat the products separately, this is the
    // "balanced" equation, where substrates have negative
    // multiplicity and products have positive multiplicity.
    //
    // This allows a better treatment of "no change" species
    // such as happen in an enzymatic reaction.  Note that
    // destroying a member of a substrate species, followed by
    // creating a member of the same species, would not necessarily
    // bring you back to the same place vis a vis the parameter
    // distribution.
    std::map<species*, int> reactionDeltas;

    // This is just the sum of the substrate multiplicities, but we
    // don't want to compute that repeatedly.
    int theArity;

    // The molar reaction rate.
    double rate;
  
  public:
    /*! \brief Keeps running count of all reaction events.

    This is a dumpable quantity. */
    static int reactionEventCount;

    /*! \brief The number of reactions that exist.

    This is a dumpable quantity. */
    static int reactionCount;

    /*! \brief Keeps a running count of reaction activations.

    An activation is calculating propensity and rescheduling.  When
    a reaction happens, many reactions may need to reschedule, and these
    reschedulings are the lengthiest part of doing the reaction.

    I suspect that the rate of activation by clock time will be more
    or less constant from simulation to simulation, unlike the
    reaction rate by clock time.  If that is so, then this dumpable
    would be useful as a way of seeing how hard the simulator is
    working. */
    static double activationCount;

    // Experimental scheduling approximation/optimization.
    double lastPropensity;

    static void
    setTolerance(double tolerance)
    {
      lowSensitive = 1.0 - tolerance;
      highSensitive = 1.0 + tolerance;
    }

    // The reaction count is no longer incremented in the reaction
    // constructor; rather, it is incremented when the reaction
    // occurs the first time.
    reaction(void) :
      theArity(0),
      rate(0.0),
      lastPropensity(-1.0)
    {
      reactionCount++;
    }

    ~reaction(void)
    {}
  
    void addSubstrate(species* pSpecies,
		      int multiplicity);

    int getSubstrateCount(void);

    void addProduct (species* pSpecies,
		     int multiplicity);

    void setRate(double newRate)
    {
      rate = newRate;
    }

    /*! \brief Schedules the reaction based on its current propensity. */
    // Reaction propensity calculation requires volume, hence mzrUnit.
    void
    schedule(eventQueue& rQueue,
	     mzrUnit& rMzrUnit);

    /*! \brief What a reaction event does, given what the reaction does. */
    eventResult
    doEvent(moleculizer& rMolzer,
	    mzrUnit& refMzrUnit);

    /*! \brief Add this reaction to the sensitivity vector of each
      substrate. */
    void
    sensitizeToSubstrates(mzrUnit& rMzrUnit);

    /*! \brief Propensity function.

    Propensity is used to caculate the reaction's next likely time of
    occurrance when the population of one of the reaction's substrate
    species changes. */
    double
    propensity(mzrUnit& rMzrUnit) const;

    // Reactions appear in more than one context, explicit-reactions and
    // generated-reactions.  Also, this routine doesn't do anything
    // if this reaction hasn't happened for the first time yet (that's why
    // return void here, instead of pointer to the inserted element.)
    void
    insertElt(xmlpp::Element* pParentElt) const throw(std::exception);
  };


  /*! \ingroup reactionGroup
    \brief Function class for rescheduling reactions en masse.

    This is for rescheduling those reactions whose substrates were
    affected by an event.  This is used in reaction::doEvent,
    createEvent::doEvent, volumeEvent::doEvent, etc. */
  class scheduleReaction :
    public std::unary_function<reaction*, void>
  {
    eventQueue& rQueue;
    mzrUnit& rMzrUnit;
    
  public:
    scheduleReaction(eventQueue& rEventQueue,
		     mzrUnit& refMzrUnit) :
      rQueue(rEventQueue),
      rMzrUnit(refMzrUnit)
    {}
  
    void operator()(reaction* pReaction) const
    {
      pReaction->schedule(rQueue,
			  rMzrUnit);
    }
  };
}

#endif
