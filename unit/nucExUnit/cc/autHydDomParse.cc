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
#include "stoch/stochUnit.hh"
#include "gpa/gpaEltName.hh"
#include "nucEx/autHydFam.hh"
#include "nucEx/nucExUnitParse.hh"
#include "nucEx/nucExEltName.hh"

namespace nucEx
{
  // Parses an auto-hydrolysis reaction generator specification.
  void
  addAutHydRxnFam::operator()(xmlpp::Node* pAutHydGenNode)
    throw(std::exception)
  {
    xmlpp::Element* pAutHydGenElt
      = domUtils::mustBeElementPtr(pAutHydGenNode);

    // Extract the nucleotide-binding mod-mol
    xmlpp::Element* pModMolRefElt
      = domUtils::mustGetUniqueChild(pAutHydGenElt,
				     gpa::eltName::modMolRef);
    std::string molName
      = domUtils::mustGetAttrString(pModMolRefElt,
				    gpa::eltName::modMolRef_nameAttr);
    bnd::mol* pMol = rMolUnit.findMol(molName);
    if(0 == pMol) throw bnd::unknownMolXcpt(pModMolRefElt,
				       molName);
    bnd::modMol* pModMol = dynamic_cast<bnd::modMol*>(pMol);
    if(0 == pModMol) throw bnd::badModMolCastXcpt(pModMolRefElt,
					     pMol);

    // Extract the nucleotide-binding modificatiion site.
    xmlpp::Element* pModSiteRefElt
      = domUtils::mustGetUniqueChild(pModMolRefElt,
				     bnd::eltName::modSiteRef);
    std::string modSiteName
      = domUtils::mustGetAttrString(pModSiteRefElt,
				    bnd::eltName::modSiteRef_nameAttr);
    int modSiteNdx
      = pModMol->getModSiteNdx(modSiteName);
    if(modSiteNdx < 0) throw bnd::unknownModSiteXcpt(pModSiteRefElt,
						modSiteName,
						pModMol);

    // Get the modification representing binding to the unhydrolysed nucleotide
    // (GTP).
    xmlpp::Element* pUnhydModRefElt
      = domUtils::mustGetUniqueChild
      (pAutHydGenElt,
       eltName::unhydrolysedBoundModRef);

    std::string unhydModName
      = domUtils::mustGetAttrString
      (pUnhydModRefElt,
       eltName::unhydrolysedBoundModRef_nameAttr);

    const bnd::modification* pUnhydMod
      = rMolUnit.mustGetMod(pUnhydModRefElt,
			    unhydModName);

    // Get the modification representing binding to the hydrolysed nucleotide
    // (GDP).
    xmlpp::Element* pHydModRefElt
      = domUtils::mustGetUniqueChild(pAutHydGenElt,
				     eltName::hydrolysedBoundModRef);
    std::string hydModName
      = domUtils::mustGetAttrString
      (pHydModRefElt,
       eltName::hydrolysedBoundModRef_nameAttr);

    const bnd::modification* pHydMod
      = rMolUnit.mustGetMod(pHydModRefElt,
			    hydModName);

    // Get the phosphate species.
    xmlpp::Element* pPhosphateSpeciesRefElt
      = domUtils::mustGetUniqueChild(pAutHydGenElt,
				     gpa::eltName::phosphateSpeciesRef);
    std::string phosphateSpeciesName
      = domUtils::mustGetAttrString
      (pPhosphateSpeciesRefElt,
       gpa::eltName::phosphateSpeciesRef_nameAttr);

    stoch::stochSpecies* pPhosphate
      = rStochUnit.mustGetStochSpecies(pPhosphateSpeciesRefElt,
				       phosphateSpeciesName);

    // Get the (unary) reaction rate.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pAutHydGenElt,
				     mzr::eltName::rate);
    double rate
      = domUtils::mustGetAttrDouble(pRateElt,
				    mzr::eltName::rate_valueAttr);

    // Construct a reaction rate extrapolator.  Which mode of extrapolation to
    // use is user-selectable, but for unary reactions, no
    // extrapolation is done at all at present..  Reaction rate extrapolators
    // are memory-managed by the reaction generators to which they are
    // attached.
    autHydNoExtrap* pExtrap
      = new autHydNoExtrap(rate);

    // Construct the reaction family.
    autHydFam* pReactionFamily
      = new autHydFam(pModMol,
		      modSiteNdx,
		      pUnhydMod,
		      pHydMod,
		      pPhosphate,
		      rMzrUnit,
		      pExtrap);

    rMzrUnit.addReactionFamily(pReactionFamily);
    
    // Connect the reaction family's reaction generator to the mol feature
    // of the auto-hydrolysing mol.
    pModMol->addRxnGen(pReactionFamily->getRxnGen());
  }
}
