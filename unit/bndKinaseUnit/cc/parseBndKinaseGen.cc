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
#include "mol/smallMol.hh"
#include "mol/molXcpt.hh"
#include "plex/parseOmniPlex.hh"
#include "bndKinase/bndKinaseEltName.hh"
#include "bndKinase/bndKinaseXcpt.hh"
#include "bndKinase/bndKinaseExtrap.hh"
#include "bndKinase/bndKinaseFam.hh"
#include "bndKinase/parseBndKinaseGen.hh"

namespace bndKinase
{
  void
  parseBndKinaseGen::
  operator()(xmlpp::Node* pBndKinaseGenNode) const
    throw(mzr::mzrXcpt)
  {
    xmlpp::Element* pBndKinaseGenElt
      = domUtils::mustBeElementPtr(pBndKinaseGenNode);

    // Find the omniPlex.
    plx::parserPlex parsedPlex;
    plx::omniPlex* pOmni
      = plx::findOmni(pBndKinaseGenNode,
		      rMolUnit,
		      rPlexUnit,
		      parsedPlex);

    // Parse the reference to the substrate mod-mol instance.
    xmlpp::Element* pSubstrateModMolInstanceRefElt
      = domUtils::mustGetUniqueChild(pBndKinaseGenElt,
				     eltName::substrateModMolInstanceRef);
    std::string substrateName
      = domUtils::mustGetAttrString(pSubstrateModMolInstanceRefElt,
				    eltName::substrateModMolInstanceRef_nameAttr);

    // Get the substrate mod-mol instance in the parsed plex.
    bnd::modMol* pSubstrateModMol = dynamic_cast<bnd::modMol*>
      (parsedPlex.mustGetMolByName(pSubstrateModMolInstanceRefElt,
				   substrateName));

    // Get the index of the substrate mol in the enabling complex.
    plx::plexMolSpec substrateMolNdx
      = parsedPlex.mustGetMolNdxByName(pSubstrateModMolInstanceRefElt,
				       substrateName);

    // Parse the modification site at which the substrate mol should
    // be phosphorylated.
    xmlpp::Element* pPhosSiteRefElt
      = domUtils::mustGetUniqueChild(pSubstrateModMolInstanceRefElt,
				     eltName::phosSiteRef);
    std::string phosSiteName
      = domUtils::mustGetAttrString(pPhosSiteRefElt,
				    eltName::phosSiteRef_nameAttr);

    // Get the index of the phosphorylation site.
    //
    // This "mustGetModSiteNdx" is a very funky routine that needs
    // some attention.  Note it accepts "this" (more or less) as an
    // argument.
    int phosSiteNdx 
      = pSubstrateModMol->mustGetModSiteNdx(pPhosSiteRefElt,
					    phosSiteName,
					    pSubstrateModMol);

    // Parse the name of the phosphorylation modification.
    xmlpp::Element* pPhosModRefElt
      = domUtils::mustGetUniqueChild(pBndKinaseGenElt,
				     eltName::phosphorylatedModRef);
    std::string phosphorylatedModName
      = domUtils::mustGetAttrString(pPhosModRefElt,
				    eltName::phosphorylatedModRef_nameAttr);
    const bnd::modification* pPhosphorylatedMod
      = rMolUnit.mustGetMod(phosphorylatedModName);

    // Parse the instance name of the ATP small-mol instance.
    xmlpp::Element* pAtpInstanceRefElt
      = domUtils::mustGetUniqueChild(pBndKinaseGenElt,
				     eltName::atpSmallMolInstanceRef);
    std::string atpInstanceName
      = domUtils::mustGetAttrString
      (pAtpInstanceRefElt,
       eltName::atpSmallMolInstanceRef_nameAttr);

    // Verify that the ATP mol is a small-mol.  In particular, this
    // will guarantee that ATP has one and only one binding site.
    bnd::smallMol* pAtpSmallMol = dynamic_cast<bnd::smallMol*>
      (parsedPlex.mustGetMolByName(pAtpInstanceRefElt,
				   atpInstanceName));
    if(! pAtpSmallMol) 
      throw badSmallMolInstanceXcpt(pAtpInstanceRefElt,
				    atpInstanceName);

    // Convert the instance name of the ATP small-mol into an index.
    int atpMolNdx
      = parsedPlex.mustGetMolNdxByName(pAtpInstanceRefElt,
				       atpInstanceName);

    // Parse the name of the ADP small-mol.
    xmlpp::Element* pAdpSmallMolRef
      = domUtils::mustGetUniqueChild(pBndKinaseGenElt,
				     eltName::adpSmallMolRef);
    std::string adpSmallMolName
      = domUtils::mustGetAttrString(pAdpSmallMolRef,
				    eltName::adpSmallMolRef_nameAttr);

    // Verify that ADP really is a small-mol.  In particular, this
    // will guarantee that ADP has one and only one binding site.
    bnd::smallMol* pAdpSmallMol = dynamic_cast<bnd::smallMol*>
      (rMolUnit.mustFindMol(pAdpSmallMolRef,
			    adpSmallMolName));
    if(! pAdpSmallMol)
      throw badSmallMolXcpt(pAdpSmallMolRef,
			    adpSmallMolName);

    // Parse the rate.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pBndKinaseGenElt,
				     eltName::rate);
    double rate
      = domUtils::mustGetAttrDouble(pRateElt,
				    eltName::rate_valueAttr);

    // Construct the unary reaction rate extrapolator. (There is only
    // one choice in this moleculizer version.).
    bndKinaseNoExtrap* pExtrap = new bndKinaseNoExtrap(rate);

    // Construct the reaction family and register it for memory
    // management.
    bndKinaseFam* pFamily
      = new bndKinaseFam(rMzrUnit,
			 rPlexUnit,
			 pSubstrateModMol,
			 substrateMolNdx,
			 phosSiteNdx,
			 pPhosphorylatedMod,
			 pAdpSmallMol,
			 atpMolNdx,
			 pExtrap);
    rMzrUnit.addReactionFamily(pFamily);

    // Connect the family's reaction generator to the omniplex's feature.
    // The omniplex feature will be connected to plexFamilies
    // that satisfy the omniplex's structural query.  The omniplex
    // feature will then notify this (and other) reaction generators
    // of species that sastify the omniplex's state query.
    pOmni->getSubPlexFeature()->addRxnGen(pFamily->getRxnGen());
  }
}
