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

#include "platform.hh"
#include "mzr/moleculizer.hh"
#include "mol/mol.hh"
#include "plex/plexUnit.hh"

namespace plx
{
  // Installs a new binding feature for the given pair of structural
  // sites. Any plex that binds the two given structural bindings
  // together will need use the bindingFeature associated to the pair
  // of structural sites to inform the family of decompositions of
  // the structural binding.
  plx::bindingFeature*
  plexUnit::addBindingFeature(bnd::mol* pLeftMol,
			       int leftMolSiteSpec,
			       bnd::mol* pRightMol,
			       int rightMolSiteSpec)
  {
    structuralSite leftSite = std::make_pair(pLeftMol,
					     leftMolSiteSpec);
    structuralSite rightSite = std::make_pair(pRightMol,
					      rightMolSiteSpec);
    structuralBinding sBinding = std::make_pair(leftSite, rightSite);

    bindingFeatureMap::iterator iEntry
      = bindingFeatures.insert(std::make_pair(sBinding,
					      plx::bindingFeature())).first;

    return &(iEntry->second);
  }

  // I expect this routine to be used in the constructor for plexFamily,
  // to construct the plexFamily's feature map for bindings.
  //
  // This routine looks for the feature "in both directions."
  plx::bindingFeature*
  plexUnit::findBindingFeature(bnd::mol* pLeftMol,
				int leftMolSiteSpec,
				bnd::mol* pRightMol,
				int rightMolSiteSpec)
  {
    structuralSite leftSite = std::make_pair(pLeftMol,
					     leftMolSiteSpec);
    structuralSite rightSite = std::make_pair(pRightMol,
					      rightMolSiteSpec);
    structuralBinding sBinding = std::make_pair(leftSite, rightSite);

    bindingFeatureMap::iterator iEntry
      = bindingFeatures.find(sBinding);

    // May have to try the sites in the opposite order.
    if(iEntry == bindingFeatures.end())
      {
	sBinding = std::make_pair(rightSite, leftSite);
	iEntry = bindingFeatures.find(sBinding);
      }

    return iEntry == bindingFeatures.end()
      ? 0
      : &(iEntry->second);
  }

  // For dumping state in XML.
  void
  plexUnit::insertStateElts(xmlpp::Element* pRootElt)
    throw(mzr::mzrXcpt)
  {
    // Get the model element.
    xmlpp::Element* pModelElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     mzr::eltName::model);

    // Ensure that the tagged-species element is here.
    xmlpp::Element* pTaggedSpeciesElt
      = domUtils::mustGetUniqueChild(pModelElt,
				     mzr::eltName::taggedSpecies);

    // Insert tagged-plex-species nodes.
    recognize.insertSpecies(pTaggedSpeciesElt,
			    rMzrUnit.getMolarFactor().getFactor());
  }
}
