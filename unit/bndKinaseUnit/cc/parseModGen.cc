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

#include "mol/modMol.hh"
#include "plex/parseOmniPlex.hh"
#include "bndKinase/bndKinaseEltName.hh"
#include "bndKinase/modExtrap.hh"
#include "bndKinase/modFam.hh"
#include "bndKinase/parseModGen.hh"

namespace bndKinase
{
  void
  parseModGen::
  operator()(xmlpp::Node* pModGenNode) const
    throw(mzr::mzrXcpt)
  {
    xmlpp::Element* pModGenElt
      = domUtils::mustBeElementPtr(pModGenNode);

    // Find the omniPlex.
    plx::parserPlex parsedPlex;
    plx::omniPlex* pOmni
      = plx::findOmni(pModGenNode,
		      rMolUnit,
		      rPlexUnit,
		      parsedPlex);

    // Parse the reference to the substrate mod-mol instance.
    xmlpp::Element* pSubstrateModMolInstanceRefElt
      = domUtils::mustGetUniqueChild(pModGenElt,
				     eltName::substrateModMolInstanceRef);

    std::string substrateInstanceName
      = domUtils::mustGetAttrString(pSubstrateModMolInstanceRefElt,
				    eltName::substrateModMolInstanceRef_nameAttr);

    // Get the substrate mod-mol.
    bnd::mol* pSubstrateMol
      = parsedPlex.mustGetMolByName(pSubstrateModMolInstanceRefElt,
				    substrateInstanceName);
    bnd::modMol* pSubstrateModMol
      = bnd::mustBeModMol(pSubstrateModMolInstanceRefElt,
			  pSubstrateMol);

    // Get the index of the substrate mol instance in the enabling
    // complex.
    plx::plexMolSpec substrateMolSpec
      = parsedPlex.mustGetMolNdxByName(pSubstrateModMolInstanceRefElt,
				       substrateInstanceName);

    // Parse the modification site at which the substrate mol should
    // be modified.
    xmlpp::Element* pModSiteRefElt
      = domUtils::mustGetUniqueChild(pSubstrateModMolInstanceRefElt,
				     eltName::modSiteRef);
    std::string modSiteName
      = domUtils::mustGetAttrString(pModSiteRefElt,
				    eltName::modSiteRef_nameAttr);

    // Get the index of the modification site.  The way that
    // "mustGetModSiteNdx takes its own "this" as an argument is peculiar,
    // due to wierd, fixable, inheritance thing.
    int modSiteNdx
      = pSubstrateModMol->mustGetModSiteNdx(pModSiteRefElt,
					    modSiteName);

    // Parse the name of the modification that is to be installed
    // at the specified modification site on the substrate mol.
    xmlpp::Element* pModRefElt
      = domUtils::mustGetUniqueChild(pModGenElt,
				     eltName::installedModRef);
    std::string installedModName
      = domUtils::mustGetAttrString(pModRefElt,
				    eltName::installedModRef_nameAttr);
    const bnd::modification* pInstalledMod
      = rMolUnit.mustGetMod(installedModName);

    // Parse the (unary) reaction rate.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pModGenElt,
				     eltName::rate);
    double rate
      = domUtils::mustGetAttrDouble(pRateElt,
				    eltName::rate_valueAttr);

    // Construct the unary reaction rate extrapolator. (Since this is
    // unary reaction, there is only once choice here.)
    modNoExtrap* pExtrap = new modNoExtrap(rate);

    // Construct the reaction family and register it for memory management.
    modFam* pFamily
      = new modFam(rMzrUnit,
		   rPlexUnit,
		   pSubstrateModMol,
		   substrateMolSpec,
		   modSiteNdx,
		   pInstalledMod,
		   pExtrap);

    // Connect the family's reaction generator to the omniplex's feature.
    // The omniplex feature will be connected to plexFamilies
    // that satisfy the omniplex's structural query.  The omniplex
    // feature will then notify this (and other) reaction generators
    // of species that sastify the omniplex's state query.
    pOmni->getSubPlexFeature()->addRxnGen(pFamily->getRxnGen());
  }
}
