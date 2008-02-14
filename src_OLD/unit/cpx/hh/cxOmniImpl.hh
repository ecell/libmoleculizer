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

#ifndef CPX_CXOMNIIMPL_H
#define CPX_CXOMNIIMPL_H

namespace cpx
{
  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  cxOmni(plexSpeciesType* pPlexSpecies,
	 const specType& rSpec) :
    fnd::featureContext<plexSpeciesType, specType>(pPlexSpecies,
						   rSpec)
  {}

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  omniPlexT*
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  getOmni(void) const
  {
    return this->getSpec().getOmni();
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  const plexIso&
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  getEmbedding(void) const
  {
    return this->getSpec().getInjection();
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  molSpec
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  translateMolSpec(molSpec specInOmni) const
  {
    return getEmbedding().forward.molMap[specInOmni];
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  bindingSpec
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  translateBindingSpec(bindingSpec specInOmni) const
  {
    return getEmbedding().forward.bindingMap[specInOmni];
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  siteSpec
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  translateSiteSpec(siteSpec specInOmni) const
  {
    return siteSpec(translateMolSpec(specInOmni.molNdx()),
			specInOmni.siteNdx());
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  plexFamilyT&
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  getPlexFamily(void) const
  {
    return this->getSpecies()->rFamily;
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  int
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  getPop(void) const
  {
    return this->getSpecies()->getPop();
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  double
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  getPlexWeight(void) const
  {
    return this->getSpecies()->getWeight();
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  const siteToShapeMap&
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  getSiteToShapeMap(void) const
  {
    return this->getSpecies()->siteToShapeMap;
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  const std::vector<molParam>&
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  getMolParams(void) const
  {
    return this->getSpecies()->molParams;
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  molParam
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  getMolParam(molSpec theMolSpec) const
  {
    const std::vector<molParam>& rMolParams = getMolParams();
    return rMolParams[theMolSpec];
  }

  template<class molT,
	   class plexSpeciesT,
	   class plexFamilyT,
	   class omniPlexT>
  const std::vector<molT*>&
  cxOmni<molT,
	 plexSpeciesT,
	 plexFamilyT,
	 omniPlexT>::
  getMols(void) const
  {
    return getPlexFamily().getParadigm().mols;
  }
}

#endif // CPX_CXOMNIIMPL_H
