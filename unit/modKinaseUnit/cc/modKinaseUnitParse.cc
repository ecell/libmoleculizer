////////////////////////////////////////////////////////////////////////////
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

#include "mzr/pchem.hh"
#include "mol/modMol.hh"
#include "mol/molUnit.hh"
#include "mol/molEltName.hh"
#include "stoch/stochSpecies.hh"
#include "stoch/stochUnit.hh"
#include "dimer/dimerEltName.hh"
#include "modKinase/nucBindFam.hh"
#include "modKinase/modKinaseFam.hh"
#include "modKinase/stochPtaseFam.hh"
#include "modKinase/modKinaseUnitParse.hh"
#include "modKinase/kinaseEltName.hh"

namespace kinase
{
  void
  addNucBindRxnFam::operator()(xmlpp::Node* pNucBindRxnGenNode) const
    throw(std::exception)
  {
    xmlpp::Element* pNucBindRxnGenElt
      = domUtils::mustBeElementPtr(pNucBindRxnGenNode);

    // Get the name of the nucleotide-binding mol.
    xmlpp::Element* pModMolRefElt
      = domUtils::mustGetUniqueChild(pNucBindRxnGenElt,
				     eltName::modMolRef);
    std::string modMolName
      = domUtils::mustGetAttrString(pModMolRefElt,
				    eltName::modMolRef_nameAttr);

    // Look up the nucleotide-binding mol, and make sure that it's
    // a modMol.
    bnd::modMol* pMol = rMolUnit.mustFindModMol(pModMolRefElt,
						modMolName);

    // Get the name of the modification site used to model
    // nucleotide binding.
    xmlpp::Element* pModSiteRefElt
      = domUtils::mustGetUniqueChild(pModMolRefElt,
				     bnd::eltName::modSiteRef);
    std::string modSiteName
      = domUtils::mustGetAttrString(pModSiteRefElt,
				    bnd::eltName::modSiteRef_nameAttr);
    int modSiteNdx
      = pMol->mustGetModSiteNdx(pModSiteRefElt,
				modSiteName,
				pMol);

    // Get the modification representing bound nucleotide at the above
    // modification site.
    xmlpp::Element* pNucBoundModRefElt
      = domUtils::mustGetUniqueChild(pNucBindRxnGenElt,
				     eltName::nucleotideBoundModRef);
    std::string nucBoundModName
      = domUtils::mustGetAttrString
      (pNucBoundModRefElt,
       eltName::nucleotideBoundModRef_nameAttr);

    // This routine should be made parallel with mustBetModMol, etc.
    const bnd::modification* pNucBoundMod
      = rMolUnit.mustGetMod(nucBoundModName);

    // Get the "none" modification.
    xmlpp::Element* pNoneModRefElt
      = domUtils::mustGetUniqueChild(pNucBindRxnGenElt,
				     eltName::noneModRef);
    std::string noneModName
      = domUtils::mustGetAttrString(pNoneModRefElt,
				    eltName::noneModRef_nameAttr);
    const bnd::modification* pNoneMod
      = rMolUnit.mustGetMod(noneModName);

    // Get the nucleotide, a stochastirator species.
    xmlpp::Element* pNucSpeciesRefElt
      = domUtils::mustGetUniqueChild(pNucBindRxnGenElt,
				     eltName::nucleotideSpeciesRef);
    std::string nucleotideSpeciesName
      = domUtils::mustGetAttrString(pNucSpeciesRefElt,
				    eltName::nucleotideSpeciesRef_nameAttr);
    stoch::stochSpecies* pNucleotide
      = rStochUnit.mustGetStochSpecies(pNucSpeciesRefElt,
				       nucleotideSpeciesName);

    // Get the on-rate.
    xmlpp::Element* pOnRateElt
      = domUtils::mustGetUniqueChild(pNucBindRxnGenElt,
				     dimer::eltName::onRate);
    double onRate
      = domUtils::mustGetAttrDouble(pOnRateElt,
				    dimer::eltName::onRate_valueAttr);

    // Get the off-rate.
    xmlpp::Element* pOffRateElt
      = domUtils::mustGetUniqueChild(pNucBindRxnGenElt,
				     dimer::eltName::offRate);
    double offRate
      = domUtils::mustGetAttrDouble(pOffRateElt,
				    dimer::eltName::offRate_valueAttr);

    // Get the rate extrapolator optional attribute, if provided.
    xmlpp::Attribute* pRateExtrapolatorAttr
      = pNucBindRxnGenElt->get_attribute
      (eltName::nucleotideBindGen_rateExtrapAttr);

    // Construct the binding rate extrapolator according to the option,
    // using the mass extrapolator as the default.
    nucBindExtrapolator* pNucBindExtrap = 0;
    if((0 != pRateExtrapolatorAttr)
       && (pRateExtrapolatorAttr->get_value()
	   == eltName::nucleotideBindGen_rateExtrap_none))
      {
	pNucBindExtrap = new nucBindNoExtrap(onRate);
      }
    else
      {
	pNucBindExtrap
	  = new nucBindMassExtrap(onRate,
				  pMol->getDefaultState()->getMolWeight(),
				  pNucleotide->getWeight());
      }

    // There is only one option for nucleotide unbinding, as usual for
    // unary reactions.
    nucUnbindExtrapolator* pNucUnbindExtrap = new nucUnbindNoExtrap(offRate);

    // Construct the reaction family.
    nucBindFam* pFamily = new nucBindFam(pMol,
					 modSiteNdx,
					 pNucBoundMod,
					 pNoneMod,
					 pNucleotide,
					 rMzrUnit,
					 pNucBindExtrap,
					 pNucUnbindExtrap);
    rMzrUnit.addReactionFamily(pFamily);

    // Attach the family's reaction generator to the mol feature
    // of the nucleotide binding mol.
    pMol->addRxnGen(pFamily->getRxnGen());
  }

  void
  setModSiteActive::operator()(xmlpp::Node* pPhosModSiteRefNode)
    throw(std::exception)
  {
    xmlpp::Element* pPhosModSiteRefElt
      = domUtils::mustBeElementPtr(pPhosModSiteRefNode);

    // This only works in more than one context because the attribute
    // is always named "name."
    std::string phosModSiteName
      = domUtils::mustGetAttrString(pPhosModSiteRefElt,
				    eltName::phosModSiteRef_nameAttr);

    int modSiteNdx = pModMol->mustGetModSiteNdx(pPhosModSiteRefElt,
						phosModSiteName,
						pModMol);
    rMask[modSiteNdx] = true;
  }

  void
  addKinaseRxnFam::operator()(xmlpp::Node* pKinaseGenNode) const
    throw(std::exception)
  {
    xmlpp::Element* pKinaseGenElt
      = domUtils::mustBeElementPtr(pKinaseGenNode);

    // Get the kinase modMol.
    xmlpp::Element* pKinaseModMolRefElt
      = domUtils::mustGetUniqueChild(pKinaseGenElt,
				     eltName::kinaseModMolRef);
    std::string kinaseModMolName
      = domUtils::mustGetAttrString(pKinaseModMolRefElt,
				    eltName::kinaseModMolRef_nameAttr);
    bnd::modMol* pKinaseModMol
      = rMolUnit.mustFindModMol(pKinaseModMolRefElt,
				kinaseModMolName);

    // Get the modification site on the kinase mod mol that is being
    // used to model ATP/ADP binding.
    xmlpp::Element* pAtpModSiteRefElt
      = domUtils::mustGetUniqueChild(pKinaseModMolRefElt,
				     eltName::atpModSiteRef);
    std::string atpModSiteName
      = domUtils::mustGetAttrString(pAtpModSiteRefElt,
				    eltName::atpModSiteRef_nameAttr);
    int atpModSiteNdx
      = pKinaseModMol->mustGetModSiteNdx(pAtpModSiteRefElt,
					 atpModSiteName,
					 pKinaseModMol);

    // Get the substrate modMol.
    xmlpp::Element* pSubstrateModMolRefElt
      = domUtils::mustGetUniqueChild(pKinaseGenElt,
				     eltName::substrateModMolRef);
    std::string substrateModMolName
      = domUtils::mustGetAttrString(pSubstrateModMolRefElt,
				    eltName::substrateModMolRef_nameAttr);
    bnd::modMol* pSubstrateModMol
      = rMolUnit.mustFindModMol(pSubstrateModMolRefElt,
				substrateModMolName);

    // Process the list of modification sites that the kinase can phosphorylate.
    xmlpp::Node::NodeList phosModSiteRefNodes
      = pSubstrateModMolRefElt->get_children(eltName::phosModSiteRef);

    std::vector<bool> activityMask(pSubstrateModMol->modSiteCount(),
				   false);

    // Set the entries corresponding to named modification sites to "true,"
    // indicating that they can be phosphorylated by the kinase.
    std::for_each(phosModSiteRefNodes.begin(),
		  phosModSiteRefNodes.end(),
		  setModSiteActive(activityMask,
				   pSubstrateModMol));

    // Get the modification that is being used to model phosphorylation.
    xmlpp::Element* pPhosphorylatedModRefElt
      = domUtils::mustGetUniqueChild(pKinaseGenElt,
				     eltName::phosphorylatedModRef);
    std::string phosphorylatedModName
      = domUtils::mustGetAttrString(pPhosphorylatedModRefElt,
				    eltName::phosphorylatedModRef_nameAttr);
    const bnd::modification* pPhosphorylatedMod
      = rMolUnit.mustGetMod(phosphorylatedModName);

    // Get the modification that is being used to model ATP binding.
    xmlpp::Element* pAtpBoundModRefElt
      = domUtils::mustGetUniqueChild(pKinaseGenElt,
				     eltName::atpBoundModRef);
    std::string atpBoundModName
      = domUtils::mustGetAttrString(pAtpBoundModRefElt,
				    eltName::atpBoundModRef_nameAttr);
    const bnd::modification* pAtpBoundMod
      = rMolUnit.mustGetMod(atpBoundModName);

    // Get the modification that is being used to model ADP binding.
    xmlpp::Element* pAdpBoundModRefElt
      = domUtils::mustGetUniqueChild(pKinaseGenElt,
				     eltName::adpBoundModRef);
    std::string adpBoundModName
      = domUtils::mustGetAttrString(pAdpBoundModRefElt,
				    eltName::adpBoundModRef_nameAttr);
    const bnd::modification* pAdpBoundMod
      = rMolUnit.mustGetMod(adpBoundModName);

    // Get the unmodification.
    xmlpp::Element* pNoneModRefElt
      = domUtils::mustGetUniqueChild(pKinaseGenElt,
				     eltName::noneModRef);
    std::string noneModName
      = domUtils::mustGetAttrString(pNoneModRefElt,
				    eltName::noneModRef_nameAttr);
    const bnd::modification* pNoneMod
      = rMolUnit.mustGetMod(noneModName);

    // Get the (binary) reaction rate.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pKinaseGenElt,
				     mzr::eltName::rate);
    double rate
      = domUtils::mustGetAttrDouble(pRateElt,
				    mzr::eltName::rate_valueAttr);

    // Get the optional rateExtrapolator attribute.
    xmlpp::Attribute* pRateExtrapolatorAttr
      = pKinaseGenElt->get_attribute
      (eltName::kinaseGen_rateExtrapAttr);

    // Construct reaction rate extrapolator according to the option
    // chosen.
    modKinaseExtrapolator *pExtrap = 0;
    if((0 != pRateExtrapolatorAttr)
       && (pRateExtrapolatorAttr->get_value()
	   == eltName::kinaseGen_rateExtrap_none))
      {
	pExtrap = new modKinaseNoExtrap(rate);
      }
    else
      {
	pExtrap = new modKinaseMassExtrap
	  (rate,
	   pKinaseModMol->getDefaultState()->getMolWeight(),
	   pSubstrateModMol->getDefaultState()->getMolWeight());
      }
				
    // Construct the reaction family.
    modKinaseFam* pFamily = new modKinaseFam(pKinaseModMol,
					     atpModSiteNdx,
					     pSubstrateModMol,
					     pPhosphorylatedMod,
					     pAtpBoundMod,
					     pAdpBoundMod,
					     pNoneMod,
					     activityMask,
					     rMzrUnit,
					     pExtrap);
    rMzrUnit.addReactionFamily(pFamily);

    // Connect the fammily's reaction generators to the kinase and substrate
    // mols.
    pKinaseModMol->addRxnGen(pFamily->getKinaseRxnGen());
    pSubstrateModMol->addRxnGen(pFamily->getSubstrateRxnGen());
  }

  void
  addPtaseRxnFam::operator()(xmlpp::Node* pPtaseGenNode) const
    throw(std::exception)
  {
    xmlpp::Element* pPtaseGenElt
      = domUtils::mustBeElementPtr(pPtaseGenNode);

    // Look up the substrate modMol.
    xmlpp::Element* pSubstrateModMolRefElt
      = domUtils::mustGetUniqueChild(pPtaseGenElt,
				     eltName::substrateModMolRef);
    std::string substrateModMolName
      = domUtils::mustGetAttrString(pSubstrateModMolRefElt,
				    eltName::substrateModMolRef_nameAttr);
    bnd::modMol* pSubstrateModMol
      = rMolUnit.mustFindModMol(pSubstrateModMolRefElt,
				substrateModMolName);

    // Process the list of modification sites that the phosphatase can
    // dephosphorylate.
    std::vector<bool> activityMask(pSubstrateModMol->modSiteCount(),
				   false);
    xmlpp::Node::NodeList phosModSiteRefNodes
      = pSubstrateModMolRefElt->get_children(eltName::phosModSiteRef);
    std::for_each(phosModSiteRefNodes.begin(),
		  phosModSiteRefNodes.end(),
		  setModSiteActive(activityMask,
				   pSubstrateModMol));

    // Get the phosphatase, a stochastirator species.
    xmlpp::Element* pPtaseStochSpeciesRefElt
      = domUtils::mustGetUniqueChild(pPtaseGenElt,
				     eltName::ptaseStochSpeciesRef);
    std::string ptaseStochSpeciesName
      = domUtils::mustGetAttrString(pPtaseStochSpeciesRefElt,
				    eltName::ptaseStochSpeciesRef_nameAttr);
    stoch::stochSpecies* pPtase
      = rStochUnit.mustGetStochSpecies(pPtaseStochSpeciesRefElt,
				       ptaseStochSpeciesName);

    // Get phosphate, a stochastirator species.
    xmlpp::Element* pPhosphateSpeciesRefElt
      = domUtils::mustGetUniqueChild(pPtaseGenElt,
				     eltName::phosphateSpeciesRef);
    std::string phosphateSpeciesName
      = domUtils::mustGetAttrString(pPhosphateSpeciesRefElt,
				    eltName::phosphateSpeciesRef_nameAttr);
    stoch::stochSpecies* pPhosphate
      = rStochUnit.mustGetStochSpecies(pPhosphateSpeciesRefElt,
				       phosphateSpeciesName);

    // Get the modification representing phophorylation.
    xmlpp::Element* pPhosphorylatedModRefElt
      = domUtils::mustGetUniqueChild(pPtaseGenElt,
				     eltName::phosphorylatedModRef);
    std::string phosphorylatedModName
      = domUtils::mustGetAttrString(pPhosphorylatedModRefElt,
				    eltName::phosphorylatedModRef_nameAttr);
    const bnd::modification* pPhosphorylatedMod
      = rMolUnit.mustGetMod(phosphorylatedModName);

    // Get the unmodification.
    xmlpp::Element* pNoneModRefElt
      = domUtils::mustGetUniqueChild(pPtaseGenElt,
				     eltName::noneModRef);
    std::string noneModName
      = domUtils::mustGetAttrString(pNoneModRefElt,
				    eltName::noneModRef_nameAttr);
    const bnd::modification* pNoneMod
      = rMolUnit.mustGetMod(noneModName);

    // Get the binary reaction rate.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pPtaseGenElt,
				     mzr::eltName::rate);
    double rate
      = domUtils::mustGetAttrDouble(pRateElt,
				    mzr::eltName::rate_valueAttr);

    // Look for the optional rate extrapolator element.
    xmlpp::Attribute* pRateExtrapolatorAttr
      = pPtaseGenElt->get_attribute
      (eltName::ptaseGen_rateExtrapAttr);

    // Construct reaction rate extrapolator according to the option,
    // with mass extrapolator the default.
    stochPtaseExtrapolator* pExtrap = 0;
    if((0 != pRateExtrapolatorAttr)
       && (pRateExtrapolatorAttr->get_value()
	   == eltName::ptaseGen_rateExtrap_none))
      {
	pExtrap = new stochPtaseNoExtrap(rate);
      }
    else
      {
	pExtrap = new stochPtaseMassExtrap
	  (rate,
	   pSubstrateModMol->getDefaultState()->getMolWeight(),
	   pPtase->getWeight());
      }


    // Construct the reaction family.
    stochPtaseFam* pFamily = new stochPtaseFam(pSubstrateModMol,
					       pPtase,
					       pPhosphate,
					       pPhosphorylatedMod,
					       pNoneMod,
					       activityMask,
					       rMzrUnit,
					       pExtrap);
    rMzrUnit.addReactionFamily(pFamily);

    // Connect the reaction generator to the substrate mol.
    pSubstrateModMol->addRxnGen(pFamily->getRxnGen());
  }
}
