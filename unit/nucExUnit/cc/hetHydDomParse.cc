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
#include "mol/molEltName.hh"
#include "plex/plexDomParse.hh"
#include "plex/plexFamily.hh"
#include "plex/plexEltName.hh"
#include "stoch/stochUnit.hh"
#include "gpa/gpaEltName.hh"
#include "nucEx/hetHydFam.hh"
#include "nucEx/nucExUnitParse.hh"
#include "nucEx/nucExEltName.hh"

namespace nucEx
{
  // Parses hetero-hydrolysis-gen, creating a reaction generator for
  // the form of the nucleotide (GTP) hydrolysis reaction that is mediated
  // by an enzyme (Sst2) other than the nucleotide-binding protein (Gpa1).
  void
  addHetHydRxnFam::operator()(xmlpp::Node* pHetHydGenNode)
    throw(std::exception)
  {
    xmlpp::Element* pHetHydGenElt
      = domUtils::mustBeElementPtr(pHetHydGenNode);

    // Parse the enabling omniplex, the dimer of the helper enzyme (Sst2)
    // and the nucleotide-binding protein (Gpa1).
    xmlpp::Element* pEnablingPlexElt
      = domUtils::mustGetUniqueChild(pHetHydGenElt,
				     plx::eltName::plex);
    plx::parserPlex enablingPlex;
    plx::plexFamily* pEnablingPlexFamily
      = recognizePlexElt(pEnablingPlexElt,
			 enablingPlex,
			 rMolUnit,
			 rPlexUnit);

    // Parse the instance name of the nucleotide-binding protein in
    // the enabling complex.
    xmlpp::Element* pTargetModMolInstanceRefElt
      = domUtils::mustGetUniqueChild(pHetHydGenElt,
				     gpa::eltName::targetModMolInstanceRef);
    std::string targetModMolInstanceName
      = domUtils::mustGetAttrString
      (pTargetModMolInstanceRefElt,
       gpa::eltName::targetModMolInstanceRef_nameAttr);

    plx::plexMolSpec targetMolSpec
      = enablingPlex.mustGetMolNdxByName(pTargetModMolInstanceRefElt,
					 targetModMolInstanceName);

    // Perhaps something like "mustBeModMolInstance" in the plexUnit
    // would be good, giving better diagnostic.
    bnd::modMol* pTargetModMol
      = bnd::mustBeModMolPtr(pTargetModMolInstanceRefElt,
			     enablingPlex.mols[targetMolSpec]);

    // Get the name of the modification site being used to model
    // nucleotide binding.
    xmlpp::Element* pModSiteRefElt
      = domUtils::mustGetUniqueChild(pTargetModMolInstanceRefElt,
				     bnd::eltName::modSiteRef);
    std::string modSiteName
      = domUtils::mustGetAttrString(pModSiteRefElt,
				    bnd::eltName::modSiteRef_nameAttr);
    int modSiteNdx
      = pTargetModMol->mustGetModSiteNdx(pModSiteRefElt,
					 modSiteName,
					 pTargetModMol);

    // Get the modification that represents the bound, unhydrolysed
    // nucleotide (GTP).
    xmlpp::Element* pUnhydrolysedBoundModRefElt
      = domUtils::mustGetUniqueChild
      (pHetHydGenElt,
       eltName::unhydrolysedBoundModRef);

    std::string unhydrolysedBoundModName
      = domUtils::mustGetAttrString
      (pUnhydrolysedBoundModRefElt,
       eltName::unhydrolysedBoundModRef_nameAttr);

    const bnd::modification* pUnhydrolysedBoundMod
      = rMolUnit.mustGetMod(pUnhydrolysedBoundModRefElt,
			    unhydrolysedBoundModName);

    // Get the modification that represents the bound, hydrolysed
    // nucleotide (GDP).
    xmlpp::Element* pHydrolysedBoundModRefElt
      = domUtils::mustGetUniqueChild(pHetHydGenElt,
				     eltName::hydrolysedBoundModRef);
    std::string hydrolysedBoundModName
      = domUtils::mustGetAttrString
      (pHydrolysedBoundModRefElt,
       eltName::hydrolysedBoundModRef_nameAttr);

    const bnd::modification* pHydrolysedBoundMod
      = rMolUnit.mustGetMod(pHydrolysedBoundModRefElt,
			    hydrolysedBoundModName);

    // Get the stochSpecies phosphate, which is a "by product."
    xmlpp::Element* pPhosphateSpeciesRefElt
      = domUtils::mustGetUniqueChild(pHetHydGenElt,
				     gpa::eltName::phosphateSpeciesRef);
    std::string phosphateSpeciesName
      = domUtils::mustGetAttrString
      (pPhosphateSpeciesRefElt,
       gpa::eltName::phosphateSpeciesRef_nameAttr);

    stoch::stochSpecies* pPhosphateSpecies
      = rStochUnit.mustGetStochSpecies(pPhosphateSpeciesRefElt,
				       phosphateSpeciesName);

    // Get the (unary) reaction rate.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pHetHydGenElt,
				     mzr::eltName::rate);
    double hydrolysisRate
      = domUtils::mustGetAttrDouble(pRateElt,
				    mzr::eltName::rate_valueAttr);

    // Construct a reaction rate extrapolator.  Which reaction rate
    // extrapolation method to use is user-selectable, but there is only
    // one rate extrapolator at present for unary reactions.
    // Reaction rate extrapolators are memory-managed by the reaction
    // generators to which they are attached.
    hetHydExtrapolator* pExtrap
      = new hetHydNoExtrap(hydrolysisRate);

    // Construct the reaction family.
    hetHydFam* pFamily = new hetHydFam(pTargetModMol,
				       targetMolSpec,
				       modSiteNdx,
				       pUnhydrolysedBoundMod,
				       pHydrolysedBoundMod,
				       pPhosphateSpecies,
				       rMzrUnit,
				       pExtrap);
    rMzrUnit.addReactionFamily(pFamily);
    

    // Connect the reaction family's reaction generator to the
    // omniplex feature of the enabling complex.
    pEnablingPlexFamily->getSubPlexFeature()->addRxnGen(pFamily->getRxnGen());
  }
}
