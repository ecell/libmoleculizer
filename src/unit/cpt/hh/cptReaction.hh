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

#ifndef CPTREACTION_H
#define CPTREACTION_H

#include "fnd/event.hh"
#include "fnd/basicReaction.hh"
#include "cpt/compartmentSpecies.hh"
#include "fnd/sensitive.hh"
#include "cpt/propensityDistro.hh"

namespace cpt
{
  class globalReaction;
  class cptApp;
  class propensityDistro;

  // In this case, the "stimulus" argument to "respond" doesn't say more about
  // the character of the event, just packages up auxiliary arguments.
  //
  // It might also be necessary to package up auxiliary arguments like this in
  // cases where there really is a "stimulus" part in the argument.
  //
  // (Also working here is the assumption in the definition of "sensitive"
  // that the "stimulus" argument is const; we can't just pass cptApp&
  // as the "stimulus."
  class reactionStim
  {
  public:
    propensityDistro& rDistro;

    reactionStim(propensityDistro& refDistro) :
      rDistro(refDistro)
    {}
  };
  
  // This class of reactions assumes that all of its reactants are
  // compartmentSpecies in a single compartment.
  //
  // A globalReaction has one of these in every compartment.
  class cptReaction :
    public fnd::basicReaction<compartmentSpecies>,
    public fnd::event<cptApp>,
    public fnd::sensitive<reactionStim>
  {
    globalReaction& rGlobalReaction;
    int compartmentIndex;

    // Reactions remember their propensity from one calculation to the next,
    // so that, when this reaction's propensity is recalculated, updating the
    // total propensity of all reactions can be a local operation.
    //    double lastCalculatedPropensity;

    // Reaction needs to remember where it is in the propensity distribution
    // so that it can be moved to the front of the distribution when the
    // reaction happens.  This is intended to be a "lazy" way of keeping the
    // more probable reactions toward the front of the distribution.
    propensityDistro::iterator iDistroEntry;

    static int
    eventCount;
    
  public :
    cptReaction(globalReaction& refGlobalReaction,
		int compartmentNdx);

    globalReaction&
    getGlobalReaction(void) const
    {
      return rGlobalReaction;
    }

    int
    getCompartmentIndex(void) const
    {
      return compartmentIndex;
    }

    compartment*
    getCompartment(void) const;

    static int
    getEventCount(void)
    {
      return eventCount;
    }

    // This also adds the reaction to the sensitivity list of the added
    // reactant.
    void
    addReactant(compartmentSpecies* pSpecies,
		int multiplicity);

    // Returns the current propensity.
    double
    propensity(void) const;

    void
    setDistroEntry(propensityDistro::iterator iEntry)
    {
      iDistroEntry = iEntry;
    }

    fnd::eventResult
    happen(cptApp& rApp)
      throw(std::exception);

    // Updates the reaction's entry in the propensity distro, as well as the
    // propensity distro's total propensity.
    void
    respond(const reactionStim& rStim);
  };
}

#endif // CPTREACTION_H
