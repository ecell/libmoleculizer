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

#ifndef STATEVAR_H
#define STATEVAR_H

/*! \file stateVar.hh
  \ingroup speciesGroup
  \brief Defines base class for moleculizer state variables. */

#include <vector>
#include <set>

namespace mzr
{
  class reaction;

  /*! \ingroup speciesGroup
    \brief Base class for moleculizer state variables.

    This is the base class not only for species, but also for global
    state variables like volume, and maybe someday temperature. */
  class stateVar
  {
    /*! \brief Reactions sensitive to this state variable.

      These might more appropriately be called "concurrent processes" at this
      level of generality. */
    std::vector<reaction*> sensitiveReactions;

  public:

    /*! \brief Add reactions to sensitivity list of this state variable.

    When a new reaction is created, it uses this routine to tell the
    state variable that the new reaction is listening to its state
    variable's value. */
    void
    addSensitiveReaction(reaction* pReaction)
    {
      sensitiveReactions.push_back(pReaction);
    }

    /*! \brief Insert sensitive reactions into affected reactions set.

    When a reaction updates a species, it uses this routine to keep
    track of which reactions will need to be rescheduled as a
    consequence. */
    void
    getSensitiveReactions(std::set<reaction*>& rAffectedReactions) const
    {
      rAffectedReactions.insert(sensitiveReactions.begin(),
				sensitiveReactions.end());
    }
  };
}

#endif
