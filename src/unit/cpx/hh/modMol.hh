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

#ifndef CPX_MODMOL_H
#define CPX_MODMOL_H

#include "cpx/alloMol.hh"
#include "cpx/modMolMixin.hh"
#include "cpx/modMolState.hh"
#include "cpx/stateMol.hh"

namespace cpx
{
  template<class baseMolT>
  class modMol :
    public alloMol<stateMol<baseMolT, modMolState> >, 
    public modMolMixin
  {
  public:
    typedef stateMol<baseMolT, modMolState> stateMolType;

    // Use molUnit::getModMap to convert a
    // map<string, string> into a map<string, const modification*>
    // as a preliminary to using this constructor.
    modMol(const baseMolT& rBaseMolT,
	   double molecularWeight,
	   const std::map<std::string, const modification*>& rModMap);

    // Still using the default state to get the molecular weight.
    // 
    // Use molUnit::getModMap to convert a
    // map<string, string> into a map<string, const modification*>
    // as a preliminary to using this function.
    const modMolState*
    internModMap(const std::map<std::string, const modification*>& rModMap);

    // Used to generate instance names for mod-mols in complexes in
    // state dump.
    std::string
    genInstanceName(int molInstanceNdx) const;
  };
}

#include "cpx/modMolImpl.hh"

#endif // CPX_MODMOL_H
