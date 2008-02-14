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

#ifndef CPTSTATEVAR_H
#define CPTSTATEVAR_H

#include "fnd/sensitivityList.hh"

namespace fnd
{
  template<class valueT, class eventT>
  class stateVar :
    public sensitivityList<eventT>
  {
    valueT theValue;

  public:

    typedef eventT sensitiveEventType;
    typedef valueT valueType;

    stateVar(const valueType& rInitialValue) :
      theValue(rInitialValue)
    {}

    const valueType&
    getValue(void) const
    {
      return theValue;
    }

    void
    setValue(const valueType& rNewValue)
    {
      theValue = rNewValue;
    }
    
    void
    updateValue(const valueType& rNewValue,
		sensitivityList<eventT>& rAffectedEvents)
    {
      setValue(rNewValue);
      getSensitives(rAffectedEvents);
    }
  };
}

#endif // CPTSTATEVAR_H
