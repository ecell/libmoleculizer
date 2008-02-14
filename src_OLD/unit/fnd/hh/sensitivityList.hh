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

#ifndef CPT_SENSITIVITYLIST_H
#define CPT_SENSITIVITYLIST_H

#include <vector>
#include <set>
#include <functional>
#include <algorithm>

namespace fnd
{
  // Template sensitivity list for state variables.
  // 
  // An alternative here might be just public inheritance from
  // std::set<sensitiveT*>.
  //
  // In original moleculizer, when a species's population is updated, the
  // species does not directly notify reactions having it as a substrate.  This
  // would result in each reaction's rescheduling once for each of its changed
  // substrates, instead of once when any of its substrates has changed.
  // 
  // Instead, the reaction accumulates all the reaction events that need
  // rescheduling into a set, then reschedules them all.  This avoids multiple
  // rescheduling of reactions, at the expense of forming the set.
  // 
  // In compartmental moleculizer, reactions don't reschedule, they just
  // recalculate their propensities, a less costly operation. But the same
  // technique is still used.
  template<class sensitiveT>
  class sensitivityList :
    public std::set<sensitiveT*>
  {
  public:
    typedef sensitiveT sensitiveType;

    // Here is the reason for using a set here: a reaction can add itself to a
    // reactant's sensitivity list more than once with no subsequent loss of
    // efficiency.
    void
    addSensitive(sensitiveType* pSensitive)
    {
      insert(pSensitive);
    }

    // This might be used to elicit immediate responses from the
    // sensitive items.
    template<class sensitiveFunction>
    void
    forEachSensitive(const sensitiveFunction& rFunction) const
    {
      std::for_each(this->begin(),
		    this->end(),
		    rFunction);
    }

    // For accumulating the union of several sensitivity lists
    // before eliciting the responses from the sensitive items.
    //
    // This is used to accumulate the union of the sensitivity lists
    // of the species whose populations are changed by a reaction event.
    void
    getSensitives(sensitivityList& rSensitives) const
    {
      rSensitives.insert(this->begin(),
			 this->end());
    }
  };
}

#endif // CPT_SENSITIVITYLIST_H
