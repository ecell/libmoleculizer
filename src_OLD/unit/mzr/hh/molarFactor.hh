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

#include "fnd/physConst.hh"
#include "fnd/stateVar.hh"
#include "mzr/mzrReaction.hh"

namespace mzr
{
  class molarFactorGlobal :
    public fnd::stateVar<double, mzrReaction>
  {
  public:
    // Default volume is 1 liter, so default molar factor is Avogadro's number.
    molarFactorGlobal(double initialVolume = 1.0) :
      fnd::stateVar<double, mzrReaction>(initialVolume * fnd::avogadrosNumber)
    {}

    const double
    getFactor(void) const
    {
      return getValue();
    }

    const double
    getVolume(void) const
    {
      return getFactor() / fnd::avogadrosNumber;
    }

    void
    updateFactor(double newMolarFactor,
		 fnd::sensitivityList<mzrReaction>& rAffectedReactions)
    {
      updateValue(newMolarFactor,
		  rAffectedReactions);
    }

    void
    updateVolume(double newVolume,
		 fnd::sensitivityList<mzrReaction>& rAffectedReactions)
    {
      updateValue(newVolume * fnd::avogadrosNumber,
		  rAffectedReactions);
    }
  };

  class volumeGlobal :
    public fnd::stateVar<double, mzrReaction>
  {
  public:
    volumeGlobal(double initialVolume = 1.0) :
      fnd::stateVar<double, mzrReaction>(initialVolume)
    {}

    const double
    getFactor(void) const
    {
      return getVolume() * fnd::avogadrosNumber;
    }

    const double
    getVolume(void) const
    {
      return getValue();
    }

    void
    updateFactor(double newMolarFactor,
		 fnd::sensitivityList<mzrReaction>& rAffectedReactions)
    {
      updateValue(newMolarFactor / fnd::avogadrosNumber,
		  rAffectedReactions);
    }

    void
    updateVolume(double newVolume,
		 fnd::sensitivityList<mzrReaction>& rAffectedReactions)
    {
      updateValue(newVolume,
		  rAffectedReactions);
    }
  };
}

#endif // MOLARFACTOR_H
