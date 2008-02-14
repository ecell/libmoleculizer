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

#ifndef CPX_CXMOLIMPL_H
#define CPX_CXMOLIMPL_H

namespace cpx
{
  template<class plexSpeciesT, class plexFamilyT>
  cxMol<plexSpeciesT, plexFamilyT>::
  cxMol(plexSpeciesType* pPlexSpecies,
	const molSpec& rSpec) :
    fnd::featureContext<plexSpeciesType, molSpec>(pPlexSpecies,
						      rSpec)
  {}

  template<class plexSpeciesT, class plexFamilyT>
  molSpec
  cxMol<plexSpeciesT, plexFamilyT>::
  getMolSpec(void) const
  {
    return this->getSpec();
  }

  template<class plexSpeciesT, class plexFamilyT>
  int
  cxMol<plexSpeciesT, plexFamilyT>::
  getPop(void) const
  {
    return this->getSpecies()->getPop();
  }

  template<class plexSpeciesT, class plexFamilyT>
  double
  cxMol<plexSpeciesT, plexFamilyT>::
  getPlexWeight(void) const
  {
    return this->getSpecies()->getWeight();
  }

  template<class plexSpeciesT, class plexFamilyT>
  plexFamilyT&
  cxMol<plexSpeciesT, plexFamilyT>::
  getPlexFamily(void) const
  {
    return this->getSpecies()->rFamily;
  }

  template<class plexSpeciesT, class plexFamilyT>
  const siteToShapeMap&
  cxMol<plexSpeciesT, plexFamilyT>::
  getSiteToShapeMap(void) const
  {
    return this->getSpecies()->siteToShapeMap;
  }

  template<class plexSpeciesT, class plexFamilyT>
  const std::vector<molParam>&
  cxMol<plexSpeciesT, plexFamilyT>::
  getMolParams(void) const
  {
    return this->getSpecies()->molParams;
  }

  template<class plexSpeciesT, class plexFamilyT>
  molParam
  cxMol<plexSpeciesT, plexFamilyT>::
  getMolParam(void) const
  {
    return getMolParams()[getMolSpec()];
  }
}

#endif // CPX_CXMOLIMPL_H
