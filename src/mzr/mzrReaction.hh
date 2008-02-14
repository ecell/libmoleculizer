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

#include "fnd/gillspReaction.hh"
#include "fnd/sensitive.hh"
#include "fnd/sensitivityList.hh"
#include "utl/dom.hh"
#include "mzr/mzrEvent.hh"

namespace mzr
{
  class mzrSpecies;
  class moleculizer;

  // This "message" is sent to a reaction when the reaction should
  // reschedule itself in the reaction queue.
  // 
  // An alternative to this is to make the argument in
  // sensitive<stimulusType>::respond(const stimulusType& rStimulus)
  // non-const.  Also, rMoleculizer seems "odd" as a stimulus;
  // it's really like an auxiliary parameter of a void stimulus.
  class mzrReactionStimulus
  {
  public:
    moleculizer& rMzr;

    mzrReactionStimulus(moleculizer& rMoleculizer) :
      rMzr(rMoleculizer)
    {}
  };

  class mzrReaction :
    public fnd::gillspReaction<mzrSpecies>,
    public fnd::sensitive<mzrReactionStimulus>,
    public mzrEvent
  {
    // Support for global state variables: variables to which
    // all reactions are sensitive, such as volume and temperature.
    typedef fnd::sensitivityList<mzrReaction> globalVar;
    
    class sensitizeToGlobal :
      public std::unary_function<globalVar*, void>
    {
      mzrReaction* pRxn;
    public:
      sensitizeToGlobal(mzrReaction* pReaction) :
	pRxn(pReaction)
      {}
      void
      operator()(globalVar* pGlobal) const
      {
	pGlobal->addSensitive(pRxn);
      }
    };

    // Support for "tolerance" optimization.
    double lastPropensity;
    static double lowSensitive;
    static double highSensitive;

  public:
    /*! \brief Keeps running count of all reaction events.
      This is a dumpable quantity. */
    static int reactionEventCount;

    /*! \brief The number of reactions that exist.
      This is a dumpable quantity. */
    static int reactionCount;

    // Support for tolerance optimization.
    static void
    setTolerance(double tolerance)
    {
      lowSensitive = 1.0 - tolerance;
      highSensitive = 1.0 + tolerance;
    }

    // This constructor enforces sensitization to global state variables, for
    // now just volume in Moleculizer, but e.g. temperature might also be
    // added.
    //
    // lastPropensity is set to -1 for the first time that the reaction
    // is rescheduled.
    template<class senseListIterator>
    mzrReaction(senseListIterator beginGlobalStateVars,
		senseListIterator endGlobalStateVars,
		double reactionRate = 0.0) :
      fnd::gillspReaction<mzrSpecies>(reactionRate),
      lastPropensity(-1.0)
    {
      std::for_each(beginGlobalStateVars,
		    endGlobalStateVars,
		    sensitizeToGlobal(this));

      ++reactionCount;
    }

    ~mzrReaction(void)
    {}

    // Overriding basic_reaction<mzrSpecies>::addReactant, so that
    // the sensitization happens with the right class of reaction,
    // rather than the basic_reaction template base class.
    void
    addReactant(mzrSpecies* pSpecies,
		int multiplicity);

    // Response of this reaction to message that one of its
    // reactants has changed population
    void
    respond(const mzrReactionStimulus& rStimulus);

    /*! \brief What a reaction does. */
    fnd::eventResult
    happen(moleculizer& rMolzer)
      throw(std::exception);

    // Output generation for state dump.
    xmlpp::Element*
    insertElt(xmlpp::Element* pParentElt) const 
      throw(std::exception);
  };
}

#endif
