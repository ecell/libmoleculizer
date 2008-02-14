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

#ifndef COMPARTMENTSPECIES_H
#define COMPARTMENTSPECIES_H

#include "fnd/basicSpecies.hh"

namespace cpt
{
  class globalSpecies;
  class cptReaction;
  class compartment;

  // One might not want this class to inherit from mzr::species.
  //
  // There's no reason that it should.  I expect this class to be
  // somewhat lighter-weight than a mzr::spcies.  This class is supposed to
  // contain compartment-specific information about a species, while the
  // referred to globalSpecies contains the rest.
  //
  // This class should inherit from sensitivityList, or have a
  // sensitivityList, to record sensitive reactions.
  class compartmentSpecies :
    public fnd::basicSpecies<cptReaction>
  {
    // It's not clear that this is called for.
    globalSpecies& rGlobalSpecies;
    
    // It's not clear that this is called for.
    int compartmentIndex;

  public:
    compartmentSpecies(globalSpecies& refGlobalSpecies,
		       int compartmentNdx,
		       int initialPop = 0) :
      fnd::basicSpecies<cptReaction>(initialPop),
      rGlobalSpecies(refGlobalSpecies),
      compartmentIndex(compartmentNdx)
    {}

    globalSpecies&
    getGlobalSpecies(void) const
    {
      return rGlobalSpecies;
    }
    
    int
    getCompartmentIndex(void) const
    {
      return compartmentIndex;
    }

    compartment*
    getCompartment(void) const;

    // For some reason, if I try to overload a getConc(double) member
    // function of basicSpecies with a getConc(void) here, I get "no such
    // member" error.  It doesn't like me to overload a member provided
    // by a base class.  Maybe this is as intended?
    //
    // For now, renaming the base routine getConcentration.
    double
    getConc(void) const;

    void
    notify(int notifyDepth);
  };
}

#endif // COMPARTMENTSPECIES_H
