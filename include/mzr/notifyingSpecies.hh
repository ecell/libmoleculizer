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

#ifndef NOTIFYINGSPECIES_H
#define NOTIFYINGSPECIES_H

#include "species.hh"

namespace mzr
{
  // A thing that is to be made aware when a species updates for the
  // first time.  Reaction generators are examples of these.
  template<class notifyingSpeciesType>
  class speciesNotificationTarget
  {
  public:
    virtual
    ~speciesNotificationTarget(void)
    {}

    // Typically, the notification target will construct reactions
    // that are sensitive to the notifying species.  These reactions
    // should be added to the notifying species's vector of sensitive
    // reactions here, typically using the species's addSensitiveReaction
    // method.
    virtual void
    notify(notifyingSpeciesType* pNotifier) = 0;
  };
    
  // Species that participates in the "moleculizer effect" of automatic
  // species and reaction generation.  Some do not, like stochastirator
  // species.
  class notifyingSpecies :
    public species
  {
    bool updated;
    
  public:
    notifyingSpecies(void) :
      updated(false)
    {}

    virtual
    ~notifyingSpecies(void)
    {}

    bool
    hasNotified(void) const
    {
      return updated;
    }

    // Notifies all the appropriate targets of this species first update.
    // Virtual so that each notifyingSpecies can find these targets in
    // its own way.
    virtual void
    notifyTargets(void) = 0;

    // This should be a "before" method for all inheriting classes.
    virtual void
    update(std::set<reaction*>& rAffectedReactions,
	   int delta)
    {
      // If this is the first update, notifyTargets, typically
      // constructing the list of sensitive reactions of this species.
      if(! updated)
	{
	  notifyTargets();
	  updated = true;
	}

      // Check the population and add the reactions sensitive to
      // this species to the set of affected reactions.
      species::update(rAffectedReactions,
		      delta);
    }
  };
}

#endif // NOTIFYINGSPECIES_H
