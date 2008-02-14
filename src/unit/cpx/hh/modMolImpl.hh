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

#ifndef CPX_MODMOLIMPL_H
#define CPX_MODMOLIMPL_H

namespace cpx
{
  template<class baseMolT>
  modMol<baseMolT>::
  modMol(const baseMolT& rBaseMol,
	 double molecularWeight,
	 const std::map<std::string, const modification*>& rModMap) :
    alloMol<stateMolType>(stateMolType(rBaseMol)),
    modMolMixin(rModMap)
  {
    // Intern the default state in this mol's alloMap.
    molParam defaultParam
      = this->internState(modMolState(molecularWeight,
				indexModMap(rModMap)));

    // Save pointer to the default state in the alloMap.
    this->pDefaultState = &(this->externState(defaultParam));
  }

  template<class baseMolT>
  const modMolState*
  modMol<baseMolT>::
  internModMap(const std::map<std::string, const modification*>& rModMap)
  {
    const modMolState& rDflt = *(this->getDefaultState());

    // This seems like an insane way to get the base molecular weight.
    const molState& rBaseState = rDflt;
    double baseWeight = rBaseState.getMolWeight();
    
    return this->internState(modMolState(baseWeight,
				   substituteModMap(rModMap,
						    rDflt)));
  }

  template<class baseMolT>
  std::string
  modMol<baseMolT>::
  genInstanceName(int molInstanceNdx) const
  {
    std::ostringstream oss;
    oss << "mod-mol_"
	<< molInstanceNdx;
    return oss.str();
  }
}

#endif // CPX_MODMOLIMPL_H
