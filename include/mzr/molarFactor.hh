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

#ifndef MOLARFACTOR_H
#define MOLARFACTOR_H

#include "mzr/physConst.hh"
#include "mzr/stateVar.hh"

namespace mzr
{
  // This is the state variable that handles volume.
  // 
  // In practice, the volume times Avogadro's number, or what I'm calling the
  // "molar factor" is used much more often in reaction propensity calculations.
  class molarFactor :
    public stateVar
  {
    double factor;
  public:
    // Default volume is 1 liter, so default molar factor is Avogadro's number.
    molarFactor(double volume = 1.0) :
      factor(volume * AVOGADROS_NUMBER)
    {}

    const double
    getFactor(void) const
    {
      return factor;
    }

    const double
    getVolume(void) const
    {
      return getFactor() / AVOGADROS_NUMBER;
    }

    void
    updateFactor(double newMolarFactor,
		 std::set<reaction*>& rAffectedReactions)
    {
      factor = newMolarFactor;
      stateVar::getSensitiveReactions(rAffectedReactions);
    }

    void
    updateVolume(double newVolume,
		 std::set<reaction*>& rAffectedReactions)
    {
      updateFactor(newVolume * AVOGADROS_NUMBER,
		   rAffectedReactions);
    }
  };
}

#endif // MOLARFACTOR_H
