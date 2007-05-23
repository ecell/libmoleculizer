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

#ifndef MOLSTATE_H
#define MOLSTATE_H

namespace bnd
{
    /*! \ingroup plexSpeciesGroup
    \brief (Concrete) base class for state of molecules. */
  class molState
  {
  protected:
    // Other contributions from molecular state may go into the true
    // molecular weight; hence the name.
    double baseWeight;

  public:
    molState(double molWeight) :
      baseWeight(molWeight)
    {}

    virtual ~molState(){}
  
    virtual double
    getMolWeight(void) const
    {
      return baseWeight;
    }

    bool
    operator<(const molState& rRight) const
    {
      return baseWeight < rRight.baseWeight;
    }
  };

  /*! \ingroup plexSpeciesGroup
    \brief Parameter giving the state of a mol in a complex.

    The whole reason for "interning" mol states in the mols was to save
    memory by not duplicating the state of each mol in each situation
    where the mol occurs. This strategy has the underlying limitation
    that mol states (for the time being) have to be relatively few in
    number.  It might be bad to have a mol with, say, a (really
    variable) double in its state. */
  typedef const molState* molParam;
}

#endif // MOLSTATE_H
