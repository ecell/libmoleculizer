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
#include "mzr/pchem.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/mzrUnitParse.hh"
#include "mol/molDomParse.hh"
#include "mol/molUnit.hh"
#include "plex/plexDomParse.hh"
#include "plex/plexFamily.hh"
#include "plex/plexEltName.hh"
#include "gpa/gpaEltName.hh"
#include "stoch/stochUnit.hh"
#include "nucEx/nucExFam.hh"
#include "nucEx/nucExUnitParse.hh"
#include "nucEx/nucExEltName.hh"

namespace nucEx
{
  void
  addNucExRxnFam::operator()(xmlpp::Node* pNucExRxnGenNode) throw(std::exception)
  {
    xmlpp::Element* pNucExRxnGenElt
      = domUtils::mustBeElementPtr(pNucExRxnGenNode);

    // Extract the enabling plex's element.
    //
    // A continuing thorn in the side here is that plexUnit will
    // have to parse this plex on its initial pass.  I need some way
    // for units to tell the plex unit where to look.  This could possibly
    // be done using XPATH?
    xmlpp::Element* pEnablingPlexElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     plx::eltName::plex);

    // Parse the enabling complex.
    plx::parserPlex enablingPlex;
    plx::plexFamily* pEnablingPlexFamily
      = recognizePlexElt(pEnablingPlexElt,
			 enablingPlex,
			 rMolUnit,
			 rPlexUnit);

    // Parse the instance states of the enabling complex to get the
    // precise species.  We need this to get a nominal molecular weight.
    // Note that the plexUnit ought to be done with setting up omniplexes
    // by now, so that constructing a plexSpecies should not be a problem.
    //
    // There is still the issue of units communicating to the plex unit the
    // location of plexes and omniplexes in their modeling constructs.
    // Possibly XPATHs could be used.
    //
    // This code is borrowed from plexDomParse.cc; it's a modified
    // form of "processPlexSpecies".
    std::vector<bnd::molParam> molParams
      = pEnablingPlexFamily->makeDefaultMolParams();
    xmlpp::Element* pInstanceStatesElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     plx::eltName::instanceStates);
    xmlpp::Node::NodeList modMolInstanceRefs
      = pInstanceStatesElt->get_children(plx::eltName::modMolInstanceRef);
    std::for_each(modMolInstanceRefs.begin(),
		  modMolInstanceRefs.end(),
		  plx::replaceModMolDefaultState(molParams,
						 pEnablingPlexFamily,
						 enablingPlex,
						 rMolUnit));
    // ("insertMember" here is the hack that made it possible to
    // do away with "plexDefinitions" and all that. It puts the
    // species in the family, but doesn't notify features of the
    // new species.  The plexUnit notifies features of all the
    // plexSpecies in "prepareToRun".)
    plx::plexSpecies* pEnablingSpecies
      = pEnablingPlexFamily->getMember(molParams);

    // Get the index of the mol that will undergo nucleotide exchange.
    xmlpp::Element* pTargetModMolInstanceRefElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     gpa::eltName::targetModMolInstanceRef);
    std::string modMolInstanceName
      = domUtils::mustGetAttrString
      (pTargetModMolInstanceRefElt,
       gpa::eltName::targetModMolInstanceRef_nameAttr);

    // Get the target mol as a basic mol, as well as its index
    // in the complex.
    plx::plexMolSpec targetMolSpec
      = enablingPlex.getMolNdxByName(modMolInstanceName);
    if(targetMolSpec < 0)
      throw plx::unknownMolInstanceXcpt(pTargetModMolInstanceRefElt,
				   modMolInstanceName);
    bnd::mol* pTargetMol = enablingPlex.mols[targetMolSpec];
    
    // Make sure that the target mol is a mod-mol.
    bnd::modMol* pTargetModMol = dynamic_cast<bnd::modMol*>(pTargetMol);
    if(0 == pTargetModMol)
      throw bnd::badModMolCastXcpt(pTargetModMolInstanceRefElt,
						   pTargetMol);

    // Get the name of the nucleotide binding site on the target mol.
    xmlpp::Element* pModSiteRefElt
      = domUtils::mustGetUniqueChild(pTargetModMolInstanceRefElt,
				     bnd::eltName::modSiteRef);
    std::string modSiteName
      = domUtils::mustGetAttrString(pModSiteRefElt,
				    bnd::eltName::modSiteRef_nameAttr);

    // Convert the modification site name into a modification site index.
    int nucModSiteNdx = pTargetModMol->getModSiteNdx(modSiteName);
    if(nucModSiteNdx < 0) throw bnd::unknownModSiteXcpt(pModSiteRefElt,
						   modSiteName,
						   pTargetModMol);

    // Get the modification that encodes the bound nucleotide.
    xmlpp::Element* pNucleotideBoundModRefElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     eltName::nucleotideBoundModRef);
    std::string nucleotideBoundModName
      = domUtils::mustGetAttrString
      (pNucleotideBoundModRefElt,
       eltName::nucleotideBoundModRef_nameAttr);

    const bnd::modification* pNucBoundMod
      = rMolUnit.mustGetMod(nucleotideBoundModName);

    // Get the modification that encodes no modification.
    xmlpp::Element* pNoneModRefElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     eltName::noneModRef);
    std::string noneModName
      = domUtils::mustGetAttrString(pNoneModRefElt,
				    eltName::noneModRef_nameAttr);
    const bnd::modification* pNoneMod
      = rMolUnit.mustGetMod(noneModName);

    // Get the nucleotide species.
    xmlpp::Element* pNucleotideSpeciesRefElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     eltName::nucleotideSpeciesRef);
    std::string nucSpeciesName
      = domUtils::mustGetAttrString
      (pNucleotideSpeciesRefElt,
       eltName::nucleotideSpeciesRef_nameAttr);

    stoch::stochSpecies* pStochNuc
      = rStochUnit.mustGetStochSpecies(pNucleotideSpeciesRefElt,
					      nucSpeciesName);

    // Get the enabled on-rate.
    xmlpp::Element* pEnabledOnRateElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     eltName::enabledOnRate);
    double enabledOnRate
      = domUtils::mustGetAttrDouble(pEnabledOnRateElt,
				    eltName::enabledOnRate_valueAttr);

    // Check the reaction rate extrapolation option for this reaction
    // generator.
    xmlpp::Attribute* pRateExtrapolatorAttr
      = pNucExRxnGenElt->get_attribute
      (eltName::nucleotideExchangeGen_rateExtrapAttr);

    // Construct "enabled" reaction rate extrapolator according to the option.
    // Reaction rate extrapolators are memory-managed by the reaction
    // generators that they are attached to.
    nucExBindExtrapolator* pEnabledBindExtrapolator = 0;
    if((0 != pRateExtrapolatorAttr)
       && (pRateExtrapolatorAttr->get_value()
	   == eltName::nucleotideExchangeGen_rateExtrap_none))
      {
	pEnabledBindExtrapolator
	  = new nucExBindNoExtrap(enabledOnRate);
      }
    else
      {
	pEnabledBindExtrapolator
	  = new nucExBindMassExtrap(enabledOnRate,
				    pEnablingSpecies->getWeight(),
				    pStochNuc->getWeight());
      }

    // Get the "plain" on-rate.
    xmlpp::Element* pPlainOnRateElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     eltName::plainOnRate);
    double plainOnRate
      = domUtils::mustGetAttrDouble(pPlainOnRateElt,
				    eltName::plainOnRate_valueAttr);

    
    // Construct "enabled" reaction rate extrapolator according to the option.
    // Reaction rate extrapolators are memory-managed by the reaction
    // generators that they are attached to.
    nucExBindExtrapolator* pPlainBindExtrapolator = 0;
    if((0 != pRateExtrapolatorAttr)
       && (pRateExtrapolatorAttr->get_value()
	   == eltName::nucleotideExchangeGen_rateExtrap_none))
      {
	pPlainBindExtrapolator
	  = new nucExBindNoExtrap(plainOnRate);
      }
    else
      {
	pPlainBindExtrapolator
	  = new nucExBindMassExtrap(plainOnRate,
				    pEnablingSpecies->getWeight(),
				    pStochNuc->getWeight());
      }

    // Get the enabled off-rate.  This is a unary rate, and doesn't
    // need any massaging.
    xmlpp::Element* pEnabledOffRateElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     eltName::enabledOffRate);
    double enabledOffRate
      = domUtils::mustGetAttrDouble(pEnabledOffRateElt,
				    eltName::enabledOffRate_valueAttr);

    // Construct the enabled unbinding reaction rate extrapolator.
    // At present, there is only one available rate extrapolator type for
    // this unary reaction.
    nucExUnbindExtrapolator* pEnabledUnbindExtrapolator
      = new nucExUnbindNoExtrap(enabledOffRate);

    // Get the "plain" off-rate, another unary rate.
    xmlpp::Element* pPlainOffRateElt
      = domUtils::mustGetUniqueChild(pNucExRxnGenElt,
				     eltName::plainOffRate);
    double plainOffRate
      = domUtils::mustGetAttrDouble(pPlainOffRateElt,
				    eltName::plainOffRate_valueAttr);

    // Construct the plain unbinding reaction rate extrapolator.
    // At present, there is only one available rate extrapolator type for
    // this unary reaction.
    nucExUnbindExtrapolator* pPlainUnbindExtrapolator
      = new nucExUnbindNoExtrap(plainOffRate);

    // Construct the reaction family.
    nucExFam* pReactionFamily
      = new nucExFam(pTargetModMol,
		     nucModSiteNdx,
		     pNoneMod,
		     pNucBoundMod,
		     pEnablingPlexFamily,
		     targetMolSpec,
		     pStochNuc,
		     rMzrUnit,
		     pEnabledBindExtrapolator,
		     pEnabledUnbindExtrapolator,
		     pPlainBindExtrapolator,
		     pPlainUnbindExtrapolator);
    rMzrUnit.addReactionFamily(pReactionFamily);

    // Connect the family's reaction generator to the mol feature
    // of the nucleotide exhanging mol.
    pTargetMol->addRxnGen(pReactionFamily->getRxnGen());
  }
}
