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
#include "nucEx/nucExEltName.hh"
#include "modKinase/modKinaseUnitParse.hh"
#include "modKinase/kinaseEltName.hh"
#include "scaffold/twentyElevenFam.hh"
#include "scaffold/scaffoldUnitParse.hh"
#include "scaffold/scafEltName.hh"

namespace scaf
{
  // Parses twenty-eleven-gen, creating a reaction generator for the special
  // scaffold kinase reaction in which Ste20 phosphorylates Ste11 in the
  // context of the 20:4:5:11 "enabling complex."
  void
  addTwentyElevenRxnFam::operator()(xmlpp::Node* pTwentyElevenGenNode)
    throw(std::exception)
  {
    xmlpp::Element* pTwentyElevenGenElt
      = domUtils::mustBeElementPtr(pTwentyElevenGenNode);

    // Parse the enabling omniplex.
    xmlpp::Element* pEnablingPlexElt
      = domUtils::mustGetUniqueChild(pTwentyElevenGenElt,
				     plx::eltName::plex);
    plx::parserPlex enablingPlex;
    plx::plexFamily* pEnablingPlexFamily
      = recognizePlexElt(pEnablingPlexElt,
			 enablingPlex,
			 rMolUnit,
			 rPlexUnit);

    // Parse the instance name of the kinase (Ste20) in the enabling complex.
    xmlpp::Element* pTwentyModMolInstanceRefElt
      = domUtils::mustGetUniqueChild
      (pTwentyElevenGenElt,
       eltName::twentyModMolInstanceRef);

    std::string twentyInstanceName
      = domUtils::mustGetAttrString
      (pTwentyModMolInstanceRefElt,
       eltName::twentyModMolInstanceRef_nameAttr);

    plx::plexMolSpec ste20MolSpec
      = enablingPlex.mustGetMolNdxByName(pTwentyModMolInstanceRefElt,
					 twentyInstanceName);
    bnd::modMol* pSte20
      = mustBeModMolPtr(pTwentyModMolInstanceRefElt,
			enablingPlex.mols[ste20MolSpec]);

    // Parse the name of the "atp-binding" modification site on Ste20.
    xmlpp::Element* pAtpModSiteRefElt
      = domUtils::mustGetUniqueChild(pTwentyModMolInstanceRefElt,
				     kinase::eltName::atpModSiteRef);
    std::string atpModSiteName
      = domUtils::mustGetAttrString(pAtpModSiteRefElt,
				    kinase::eltName::atpModSiteRef_nameAttr);
    int ste20AtpBindingSiteNdx
      = pSte20->mustGetModSiteNdx(pAtpModSiteRefElt,
				  atpModSiteName,
				  pSte20);

    // Parse the intance name of the "substrate" (Ste11) in the enabling
    // complex.
    xmlpp::Element* pElevenModMolInstanceRefElt
      = domUtils::mustGetUniqueChild
      (pTwentyElevenGenElt,
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

    // Parse the target phosphorylation sites on Ste11 as an "activity mask."
    xmlpp::Node::NodeList targetModSiteRefNodes
      = pElevenModMolInstanceRefElt
      ->get_children(eltName::targetModSiteRef);

    std::vector<bool> activityMask(pSte11->modSiteCount(),
				   false);

    // setModSiteActive defined and used in modKinaseUnitParse.cc
    std::for_each(targetModSiteRefNodes.begin(),
		  targetModSiteRefNodes.end(),
		  kinase::setModSiteActive(activityMask,
				   pSte11));

    // Get the modification that is being used to model ATP binding.
    xmlpp::Element* pAtpBoundModRefElt
      = domUtils::mustGetUniqueChild(pTwentyElevenGenElt,
				     kinase::eltName::atpBoundModRef);
    std::string atpBoundModName
      = domUtils::mustGetAttrString
      (pAtpBoundModRefElt,
       kinase::eltName::atpBoundModRef_nameAttr);

    const bnd::modification* pAtpBoundMod
      = rMolUnit.mustGetMod(atpBoundModName);

    // Get the modification that is being used to model ADP binding.
    xmlpp::Element* pAdpBoundModRefElt
      = domUtils::mustGetUniqueChild(pTwentyElevenGenElt,
				     kinase::eltName::adpBoundModRef);
    std::string adpBoundModName
      = domUtils::mustGetAttrString
      (pAdpBoundModRefElt,
       kinase::eltName::adpBoundModRef_nameAttr);

    const bnd::modification* pAdpBoundMod
      = rMolUnit.mustGetMod(adpBoundModName);

    // Get the modification that is being used to model phosphorylation.
    xmlpp::Element* pPhosphorylatedModRefElt
      = domUtils::mustGetUniqueChild(pTwentyElevenGenElt,
				     kinase::eltName::phosphorylatedModRef);
    std::string phosphorylatedModName
      = domUtils::mustGetAttrString
      (pPhosphorylatedModRefElt,
       kinase::eltName::phosphorylatedModRef_nameAttr);

    const bnd::modification* pPhosphorylatedMod
      = rMolUnit.mustGetMod(phosphorylatedModName);

    // Get the unmodification.
    xmlpp::Element* pNoneModRefElt
      = domUtils::mustGetUniqueChild(pTwentyElevenGenElt,
				     nucEx::eltName::noneModRef);
    std::string noneModName
      = domUtils::mustGetAttrString(pNoneModRefElt,
				    nucEx::eltName::noneModRef_nameAttr);
    const bnd::modification* pNoneMod
      = rMolUnit.mustGetMod(noneModName);

    // Get the (unary) reaction rate.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pTwentyElevenGenElt,
				     mzr::eltName::rate);
    double rate
      = domUtils::mustGetAttrDouble(pRateElt,
				    mzr::eltName::rate_valueAttr);

    // Construct a reaction rate extrapolator.  The mode of reaction rate
    // extrapolation is user-selectable, but at present there are no real
    // extrapolators for unary reactions.  Reaction rate extrapolators are
    // memory managed by the reaction rate generators to which they are
    // attached.
    twentyElevenExtrapolator* pExtrap
      = new twentyElevenNoExtrap(rate);

    // Construct the reaction family.
    twentyElevenFam* pFamily
      = new twentyElevenFam(rMzrUnit,

			    pSte20,
			    ste20AtpBindingSiteNdx,
			    ste20MolSpec,

			    pSte11,
			    activityMask,
			    ste11MolSpec,

			    pAtpBoundMod,
			    pAdpBoundMod,
			    pPhosphorylatedMod,
			    pNoneMod,

			    pExtrap);
    rMzrUnit.addReactionFamily(pFamily);

    // Connect the reaction family's reaction generator to the
    // omniplex feature of the enabling complex.
    pEnablingPlexFamily->getSubPlexFeature()->addRxnGen(pFamily->getRxnGen());
  }
}
