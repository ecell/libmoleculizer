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

#ifndef CPX_MODMOLSTATE_H
#define CPX_MODMOLSTATE_H

#include "cpx/molState.hh"
#include "cpx/modStateMixin.hh"

namespace cpx
{
  class modMolState :
    public molState, 
    public modStateMixin
  {
  public:
    modMolState(double molecularWeight,
		const modification* pModification,
		int count) :
      molState(molecularWeight),
      modStateMixin(pModification, count)
    {}

    modMolState(double molecularWeight,
		const modStateMixin& rModStateMixin) :
      molState(molecularWeight),
      modStateMixin(rModStateMixin)
    {}

    bool operator< (const modMolState& rRight) const
    {
      const molState& rThisMolState = *this;
      const molState& rRightMolState = rRight;

      if(rThisMolState < rRightMolState) return true;
      if(rRightMolState < rThisMolState) return false;

      const modStateMixin& rThisModState = *this;
      const modStateMixin& rRightModState = rRight;

      return rThisModState < rRightModState;
    }

    virtual
    ~modMolState(void)
    {}

    virtual double
    getMolWeight(void) const
    {
      return baseWeight + totalWeightDelta();
    }
  };
}

#endif // CPX_MODMOLSTATE_H
