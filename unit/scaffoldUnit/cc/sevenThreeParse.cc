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
#include "mzr/mzrEltName.hh"
#include "mol/molDomParse.hh"
#include "mol/molUnit.hh"
#include "plex/plexDomParse.hh"
#include "plex/plexFamily.hh"
#include "plex/plexEltName.hh"
#include "modKinase/modKinaseUnitParse.hh"
#include "modKinase/kinaseEltName.hh"
#include "nucEx/nucExEltName.hh"
#include "scaffold/sevenThreeFam.hh"
#include "scaffold/scaffoldUnitParse.hh"
#include "scaffold/scafEltName.hh"

namespace scaf
{
  // Parses seven-three-gen, creating reaction generator for special scaffold
  // kinase reaction in which Ste7 phosphorylated Fus3 in the context of the
  // 7:5:3 "enabling complex."
  void
  addSevenThreeRxnFam::operator()(xmlpp::Node* pSevenThreeGenNode)
    throw(std::exception)
  {
    xmlpp::Element* pSevenThreeGenElt
      = domUtils::mustBeElementPtr(pSevenThreeGenNode);

    // Parse the enabling complex.
    xmlpp::Element* pEnablingPlexElt
      = domUtils::mustGetUniqueChild(pSevenThreeGenElt,
				     plx::eltName::plex);
    plx::parserPlex enablingPlex;
    plx::plexFamily* pEnablingPlexFamily
      = recognizePlexElt(pEnablingPlexElt,
			 enablingPlex,
			 rMolUnit,
			 rPlexUnit);

    // Parse the instance name of the kinase (Ste7) in the enabling complex.
    xmlpp::Element* pSevenModMolInstanceRefElt
      = domUtils::mustGetUniqueChild(pSevenThreeGenElt,
				     eltName::sevenModMolInstanceRef);
    std::string sevenInstanceName
      = domUtils::mustGetAttrString
      (pSevenModMolInstanceRefElt,
       eltName::sevenModMolInstanceRef_nameAttr);

    plx::plexMolSpec ste7MolSpec
      = enablingPlex.mustGetMolNdxByName(pSevenModMolInstanceRefElt,
					 sevenInstanceName);
    bnd::modMol* pSte7
      = mustBeModMolPtr(pSevenModMolInstanceRefElt,
			enablingPlex.mols[ste7MolSpec]);

    // Parse the name of the "atp-binding" modification site on Ste7.
    xmlpp::Element* pAtpModSiteRefElt
      = domUtils::mustGetUniqueChild(pSevenModMolInstanceRefElt,
				     kinase::eltName::atpModSiteRef);
    std::string atpModSiteName
      = domUtils::mustGetAttrString(pAtpModSiteRefElt,
				    kinase::eltName::atpModSiteRef_nameAttr);

    // Absurd third mol argument, same as "this"?
    int ste7AtpModSiteNdx
      = pSte7->mustGetModSiteNdx(pAtpModSiteRefElt,
				 atpModSiteName,
				 pSte7);

    // Parse the name of the modification site on Ste7 whose phosphorylation is
    // a precondition of the reaction.
    xmlpp::Element* pActiveModSiteRefElt
      = domUtils::mustGetUniqueChild(pSevenModMolInstanceRefElt,
				     eltName::activeModSiteRef);
    std::string activeModSiteName
      = domUtils::mustGetAttrString(pActiveModSiteRefElt,
				    eltName::activeModSiteRef_nameAttr);
    // Absurd third argument?
    int activeModSiteNdx
      = pSte7->mustGetModSiteNdx(pActiveModSiteRefElt,
				 activeModSiteName,
				 pSte7);

    // Parse the name of the modification site on Ste7 whose phosphorylation
    // inhibits the reaction.
    xmlpp::Element* pInhibModSiteRefElt
      = domUtils::mustGetUniqueChild(pSevenModMolInstanceRefElt,
				     eltName::inhibModSiteRef);
    std::string inhibModSiteName
      = domUtils::mustGetAttrString(pInhibModSiteRefElt,
				    eltName::inhibModSiteRef_nameAttr);
    // Absurd third argument?
    int inhibModSiteNdx
      = pSte7->mustGetModSiteNdx(pInhibModSiteRefElt,
				 inhibModSiteName,
				 pSte7);

    // Parse the instance name of the "substrate" Fus3 in the enabling complex.
    xmlpp::Element* pThreeModMolInstanceRefElt
      = domUtils::mustGetUniqueChild(pSevenThreeGenElt,
				     eltName::threeModMolInstanceRef);
    std::string threeInstanceName
      = domUtils::mustGetAttrString
      (pThreeModMolInstanceRefElt,
       eltName::threeModMolInstanceRef_nameAttr);

    plx::plexMolSpec fus3MolSpec
      = enablingPlex.mustGetMolNdxByName(pThreeModMolInstanceRefElt,
					 threeInstanceName);
    bnd::modMol* pFus3
      = mustBeModMolPtr(pThreeModMolInstanceRefElt,
			enablingPlex.mols[fus3MolSpec]);

    // Parse the modification sites on Fus3 that can be phosphorylated by Ste7
    // as an "activity mask."
    xmlpp::Node::NodeList targetModSiteRefNodes
      = pThreeModMolInstanceRefElt
      ->get_children(eltName::targetModSiteRef);
  
    std::vector<bool> activityMask(pFus3->modSiteCount(),
				   false);

    // setModSiteActive defined and used in modKinaseUnitParse.cc
    std::for_each(targetModSiteRefNodes.begin(),
		  targetModSiteRefNodes.end(),
		  kinase::setModSiteActive(activityMask,
					   pFus3));

    // Get the modification that is being used to model ATP binding.
    xmlpp::Element* pAtpBoundModRefElt
      = domUtils::mustGetUniqueChild(pSevenThreeGenElt,
				     kinase::eltName::atpBoundModRef);
    std::string atpBoundModName
      = domUtils::mustGetAttrString(pAtpBoundModRefElt,
				    kinase::eltName::atpBoundModRef_nameAttr);
    const bnd::modification* pAtpBoundMod
      = rMolUnit.mustGetMod(atpBoundModName);

    // Get the modification that is being used to model ADP binding.
    xmlpp::Element* pAdpBoundModRefElt
      = domUtils::mustGetUniqueChild(pSevenThreeGenElt,
				     kinase::eltName::adpBoundModRef);
    std::string adpBoundModName
      = domUtils::mustGetAttrString(pAdpBoundModRefElt,
				    kinase::eltName::adpBoundModRef_nameAttr);
    const bnd::modification* pAdpBoundMod
      = rMolUnit.mustGetMod(adpBoundModName);

    // Get the modification that is being used to model phosphorylation.
    xmlpp::Element* pPhosphorylatedModRefElt
      = domUtils::mustGetUniqueChild(pSevenThreeGenElt,
				     kinase::eltName::phosphorylatedModRef);
    std::string phosphorylatedModName
      = domUtils::mustGetAttrString(pPhosphorylatedModRefElt,
				    kinase::eltName::phosphorylatedModRef_nameAttr);
    const bnd::modification* pPhosphorylatedMod
      = rMolUnit.mustGetMod(phosphorylatedModName);

    // Get the unmodification.
    xmlpp::Element* pNoneModRefElt
      = domUtils::mustGetUniqueChild(pSevenThreeGenElt,
				     nucEx::eltName::noneModRef);
    std::string noneModName
      = domUtils::mustGetAttrString(pNoneModRefElt,
				    nucEx::eltName::noneModRef_nameAttr);
    const bnd::modification* pNoneMod
      = rMolUnit.mustGetMod(noneModName);

    // Get the (unary) reaction rate.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pSevenThreeGenElt,
				     mzr::eltName::rate);
    double rate
      = domUtils::mustGetAttrDouble(pRateElt,
				    mzr::eltName::rate_valueAttr);

    // Construct a reaction rate extrapolator.  The mode of reaction rate
    // extrapolation is user-selectable, but at present there are no real
    // extrapolators for unary reactions.  Reaction rate extrapolators are
    // memory-managed by the reaction generators to which they are attached.
    sevenThreeExtrapolator* pExtrap
      = new sevenThreeNoExtrap(rate);

    // Construct the reaction family.
    sevenThreeFam* pFamily
      = new sevenThreeFam(rMzrUnit,

			  pSte7,
			  ste7AtpModSiteNdx,
			  ste7MolSpec,
			  activeModSiteNdx,
			  inhibModSiteNdx,

			  pFus3,
			  fus3MolSpec,
			  activityMask,

			  pAtpBoundMod,
			  pAdpBoundMod,
			  pPhosphorylatedMod,
			  pNoneMod,

			  pExtrap);

    rMzrUnit.addReactionFamily(pFamily);

    // Connect the omniplex feature of the enabling plexFamily to the reaction
    // family's reaction generator, so that the reaction generator will be
    // notified whenever a new species of complex appears that contains the
    // enabling subcomplex.
    pEnablingPlexFamily->getSubPlexFeature()->addRxnGen(pFamily->getRxnGen());
  }
}
