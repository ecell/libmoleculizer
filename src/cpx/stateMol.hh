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

#ifndef CPX_STATEMOL_H
#define CPX_STATEMOL_H

namespace cpx
{
  template<class baseMolT, class stateT>
  class stateMol :
    public baseMolT
  {
  protected:
    // The default state of the mol.  This is used as a starting point
    // for constructing other states of the mol.
    //
    // Note that this isn't set in the constructor; descendant classes
    // must set this.
    const stateT* pDefaultState;

  public:
    typedef baseMolT baseMolType;
    typedef stateT stateType;

    stateMol(const baseMolT& rBaseMol) :
      baseMolT(rBaseMol),
      pDefaultState(0)
    {}

    // Returns the default state as a reference to the actual state class.
    // One can can then use the default state as the starting point
    // for making other states of the mol.
    const stateT*
    getDefaultState(void) const
    {
      return pDefaultState;
    }

    // Note that this is/should be virtual in baseMolT, as in basicMol.
    molParam
    getDefaultParam(void) const
    {
      return pDefaultState;
    }
  };
}

#endif // CPX_STATEMOL_H
