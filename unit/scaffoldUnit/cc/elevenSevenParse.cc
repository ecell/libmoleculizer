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

#include "domUtils/domUtils.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/mzrUnitParse.hh"
#include "mol/molDomParse.hh"
#include "mol/molUnit.hh"
#include "plex/plexDomParse.hh"
#include "plex/plexFamily.hh"
#include "plex/plexEltName.hh"
#include "modKinase/kinaseEltName.hh"
#include "nucEx/nucExEltName.hh"
#include "scaffold/elevenSevenFam.hh"
#include "scaffold/scaffoldUnitParse.hh"
#include "scaffold/scafEltName.hh"

namespace scaf
{
  // Parses eleven-seven-gen, creating a reaction generator for the special
  // scaffold kinase reaction in which Ste11 phosphorylates Ste7 in the context
  // of the 11:5:7 "enabling complex."
  void
  addElevenSevenRxnFam::operator()(xmlpp::Node* pElevenSevenGenNode)
    throw(std::exception)
  {
    xmlpp::Element* pElevenSevenGenElt
      = domUtils::mustBeElementPtr(pElevenSevenGenNode);

    // Parse the enabling complex.
    xmlpp::Element* pEnablingPlexElt
      = domUtils::mustGetUniqueChild(pElevenSevenGenElt,
				     plx::eltName::plex);
    plx::parserPlex enablingPlex;
    plx::plexFamily* pEnablingPlexFamily
      = recognizePlexElt(pEnablingPlexElt,
			 enablingPlex,
			 rMolUnit,
			 rPlexUnit);

    // Parse the instance name of the kinase (Ste11) in the enabling complex.
    xmlpp::Element* pElevenModMolInstanceRefElt
      = domUtils::mustGetUniqueChild(pElevenSevenGenElt,
				     eltName::elevenModMolInstanceRef);
    std::string elevenInstanceName
      = domUtils::mustGetAttrString
      (pElevenModMolInstanceRefElt,
       eltName::elevenModMolInstanceRef_nameAttr);

    plx::plexMolSpec ste11MolSpec
      = enablingPlex.mustGetMolNdxByName(pElevenModMolInstanceRefElt,
					 elevenInstanceName);
    bnd::modMol* pSte11
      = mustBeModMolPtr(pElevenModMolInstanceRefElt,
			enablingPlex.mols[ste11MolSpec]);

    // Parse the name of the "atp-binding" modification site on Ste11.
    xmlpp::Element* pAtpModSiteRefElt
      = domUtils::mustGetUniqueChild(pElevenModMolInstanceRefElt,
				     kinase::eltName::atpModSiteRef);
    std::string atpModSiteName
      = domUtils::mustGetAttrString(pAtpModSiteRefElt,
				    kinase::eltName::atpModSiteRef_nameAttr);
    // Absurd third argument?
    int ste11AtpModSiteNdx
      = pSte11->mustGetModSiteNdx(pAtpModSiteRefElt,
				  atpModSiteName,
				  pSte11);

    // Parse the minimum phosphorylation count at which Ste11 will
    // be active.
    xmlpp::Element* pMinActivePhosCountElt
      = domUtils::mustGetUniqueChild(pElevenModMolInstanceRefElt,
				     eltName::minActivePhosCount);
    int elevenActivePhosCount
      = domUtils::mustGetAttrInt(pMinActivePhosCountElt,
				 eltName::minActivePhosCount_valueAttr);

    // Parse the "substrate" mol instance.
    xmlpp::Element* pSevenModMolInstanceRefElt
      = domUtils::mustGetUniqueChild(pElevenSevenGenElt,
				     eltName::sevenModMolInstanceRef);
    std::string sevenModMolInstanceName
      = domUtils::mustGetAttrString
      (pSevenModMolInstanceRefElt,
       eltName::sevenModMolInstanceRef_nameAttr);

    plx::plexMolSpec ste7MolSpec
      = enablingPlex.mustGetMolNdxByName(pSevenModMolInstanceRefElt,
					 sevenModMolInstanceName);
  
    bnd::modMol* pSte7
      = mustBeModMolPtr(pSevenModMolInstanceRefElt,
			enablingPlex.mols[ste7MolSpec]);

    // Parse the target modification site on Ste7, the modification site
    // that will be phosphorylated.
    xmlpp::Element* pTargetModSiteRefElt
      = domUtils::mustGetUniqueChild(pSevenModMolInstanceRefElt,
				     eltName::targetModSiteRef);
    std::string targetModSiteName
      = domUtils::mustGetAttrString(pTargetModSiteRefElt,
				    eltName::targetModSiteRef_nameAttr);
    // Absurd third argument?
    int ste7TargetModSiteNdx
      = pSte7->mustGetModSiteNdx(pTargetModSiteRefElt,
				 targetModSiteName,
				 pSte7);

    // Parse the usual modifications.
    // Get the modification that is being used to model ATP binding.
    xmlpp::Element* pAtpBoundModRefElt
      = domUtils::mustGetUniqueChild(pElevenSevenGenElt,
				     kinase::eltName::atpBoundModRef);
    std::string atpBoundModName
      = domUtils::mustGetAttrString(pAtpBoundModRefElt,
				    kinase::eltName::atpBoundModRef_nameAttr);
    const bnd::modification* pAtpBoundMod
      = rMolUnit.mustGetMod(atpBoundModName);

    // Get the modification that is being used to model ADP binding.
    xmlpp::Element* pAdpBoundModRefElt
      = domUtils::mustGetUniqueChild(pElevenSevenGenElt,
				     kinase::eltName::adpBoundModRef);
    std::string adpBoundModName
      = domUtils::mustGetAttrString(pAdpBoundModRefElt,
				    kinase::eltName::adpBoundModRef_nameAttr);
    const bnd::modification* pAdpBoundMod
      = rMolUnit.mustGetMod(adpBoundModName);

    // Get the modification that is being used to model phosphorylation.
    xmlpp::Element* pPhosphorylatedModRefElt
      = domUtils::mustGetUniqueChild(pElevenSevenGenElt,
				     kinase::eltName::phosphorylatedModRef);
    std::string phosphorylatedModName
      = domUtils::mustGetAttrString(pPhosphorylatedModRefElt,
				    kinase::eltName::phosphorylatedModRef_nameAttr);
    const bnd::modification* pPhosphorylatedMod
      = rMolUnit.mustGetMod(phosphorylatedModName);

    // Get the unmodification.
    xmlpp::Element* pNoneModRefElt
      = domUtils::mustGetUniqueChild(pElevenSevenGenElt,
				     nucEx::eltName::noneModRef);
    std::string noneModName
      = domUtils::mustGetAttrString(pNoneModRefElt,
				    nucEx::eltName::noneModRef_nameAttr);
    const bnd::modification* pNoneMod
      = rMolUnit.mustGetMod(noneModName);

    // Get the (unary) reaction rate.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pElevenSevenGenElt,
				     mzr::eltName::rate);
    double rate
      = domUtils::mustGetAttrDouble(pRateElt,
				    mzr::eltName::rate_valueAttr);

    // Construct a reaction rate extrapolator.  The mode of reaction rate
    // extrapolation is user-selectable, but at present there are no real
    // extrapolators for unary reactions.  Reaction rate extrapolators are
    // memory-managed by the reaction generators to which they are attached.
    elevenSevenExtrapolator* pExtrap
      = new elevenSevenNoExtrap(rate);

    // Construct the reaction family.
    elevenSevenFam* pFamily
      = new elevenSevenFam(rMzrUnit,

			   pSte11,
			   ste11AtpModSiteNdx,
			   ste11MolSpec,
			   elevenActivePhosCount,

			   pSte7,
			   ste7MolSpec,
			   ste7TargetModSiteNdx,

			   pAtpBoundMod,
			   pAdpBoundMod,
			   pPhosphorylatedMod,
			   pNoneMod,

			   pExtrap);

    rMzrUnit.addReactionFamily(pFamily);
    

    // Connect the family's reaction generator to the omniplex feature of the
    // enabling complex, so that the reaction generator will be notified
    // whenever a new species of complex containing the enabling subcomplex
    // appears.
    pEnablingPlexFamily->getSubPlexFeature()->addRxnGen(pFamily->getRxnGen());
  }
}
