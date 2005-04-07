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

#include "plex/parserPlex.hh"
#include "plex/parseOmniPlex.hh"
#include "mol/modMol.hh"
#include "bndKinase/bndKinaseEltName.hh"
#include "bndKinase/bndOmniGen.hh"
#include "bndKinase/bndOmniFam.hh"
#include "bndKinase/parseBndOmniGen.hh"

namespace bndKinase
{
  class parseSmallMolExchange :
    public std::unary_function<xmlpp::Node*, smallMolExchange>
  {
    bnd::molUnit& rMolUnit;
    const plx::parserPlex& rParsedPlex;
  public:
    parseSmallMolExchange(bnd::molUnit& refMolUnit,
			  const plx::parserPlex& refParsedPlex) :
      rMolUnit(refMolUnit),
      rParsedPlex(refParsedPlex)
    {}

    smallMolExchange
    operator()(xmlpp::Node* pSmallMolExchangeNode) const
      throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pSmallMolExchangeElt
	= domUtils::mustBeElementPtr(pSmallMolExchangeNode);
      
      xmlpp::Element* pSmallMolInstanceRefElt
	= domUtils::mustGetUniqueChild(pSmallMolExchangeElt,
				       eltName::smallMolInstanceRef);

      std::string smallMolInstanceName
	= domUtils::mustGetAttrString(pSmallMolInstanceRefElt,
				      eltName::smallMolInstanceRef_nameAttr);
      plx::plexMolSpec exchangedMolSpec
	= rParsedPlex.mustGetMolNdxByName(pSmallMolInstanceRefElt,
					  smallMolInstanceName);

      xmlpp::Element* pSmallMolRefElt
	= domUtils::mustGetUniqueChild(pSmallMolExchangeElt,
				       eltName::smallMolRef);

      std::string smallMolName
	= domUtils::mustGetAttrString(pSmallMolRefElt,
				      eltName::smallMolRef_nameAttr);
      
      bnd::smallMol* pReplacementMol
	= rMolUnit.mustFindSmallMol(pSmallMolInstanceRefElt,
				    smallMolName);

      return smallMolExchange(exchangedMolSpec,
			      pReplacementMol);
    }
  };

  class parseModificationExchange :
    public std::unary_function<xmlpp::Node*, modificationExchange>
  {
    bnd::molUnit& rMolUnit;
    const plx::parserPlex& rParsedPlex;
  public:
    parseModificationExchange(bnd::molUnit& refMolUnit,
			      const plx::parserPlex& refParsedPlex) :
      rMolUnit(refMolUnit),
      rParsedPlex(refParsedPlex)
    {
    }

    modificationExchange
    operator()(xmlpp::Node* pModificationExchangeNode) const
      throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pModificationExchangeElt
	= domUtils::mustBeElementPtr(pModificationExchangeNode);

      xmlpp::Element* pModMolInstanceRefElt
	= domUtils::mustGetUniqueChild(pModificationExchangeElt,
				       eltName::modMolInstanceRef);

      // Get the instance name of the mol instance in which the modification
      // exchange is to take place, and get the index of the mol in the
      // complex.
      std::string modMolInstanceName
	= domUtils::mustGetAttrString(pModMolInstanceRefElt,
				      eltName::modMolInstanceRef_nameAttr);

      plx::plexMolSpec modMolSpec
	= rParsedPlex.mustGetMolNdxByName(pModMolInstanceRefElt,
					  modMolInstanceName);

      // Get the modification site name.
      xmlpp::Element* pModSiteRefElt
	= domUtils::mustGetUniqueChild(pModMolInstanceRefElt,
				       eltName::modSiteRef);

      std::string modSiteName
	= domUtils::mustGetAttrString(pModSiteRefElt,
				      eltName::modSiteRef_nameAttr);

      // Have to use the mol to translate modification site name to
      // modification index.
      const bnd::modMol* pModMol
	= bnd::mustBeModMolPtr
	(pModMolInstanceRefElt,
	 rParsedPlex.mustGetMolByName(pModMolInstanceRefElt,
				      modMolInstanceName));
      int modSiteNdx
	= pModMol->mustGetModSiteNdx(pModSiteRefElt,
				     modSiteName);

      // Get the name of the modification that is to replace the modification
      // currently found at the modification site.
      xmlpp::Element* pInstalledModRefElt
	= domUtils::mustGetUniqueChild(pModificationExchangeElt,
				       eltName::installedModRef);
      std::string modificationName
	= domUtils::mustGetAttrString(pInstalledModRefElt,
				      eltName::installedModRef_nameAttr);

      // Look up the modification using the molUnit.
      const bnd::modification* pModification
	= rMolUnit.mustGetMod(pInstalledModRefElt,
			      modificationName);

      return modificationExchange(modMolSpec,
				  modSiteNdx,
				  pModification);
    }
  };

  void
  parseBndOmniGen::
  operator()(xmlpp::Node* pBndOmniGenNode) const
    throw(mzr::mzrXcpt)
  {
    xmlpp::Element* pBndOmniGenElt
      = domUtils::mustBeElementPtr(pBndOmniGenNode);

    // Parse the enabling omniplex.
    xmlpp::Element* pEnablingOmniElt
      = domUtils::mustGetUniqueChild(pBndOmniGenElt,
				     eltName::enablingOmniplex);
    plx::parserPlex parsedPlex;
    plx::omniPlex* pOmni
      = plx::findOmni(pEnablingOmniElt,
		      rMolUnit,
		      rPlexUnit,
		      parsedPlex);

    // Parse the small-mol exchanges.
    xmlpp::Element* pSmallMolExchangesElt
      = domUtils::mustGetUniqueChild(pBndOmniGenElt,
				     eltName::smallMolExchanges);
    xmlpp::Node::NodeList smallMolExchangeNodes
      = pSmallMolExchangesElt->get_children(eltName::smallMolExchange);
    std::vector<smallMolExchange> smallMolExchanges
      = std::vector<smallMolExchange>(smallMolExchangeNodes.size());
    std::transform(smallMolExchangeNodes.begin(),
		   smallMolExchangeNodes.end(),
		   smallMolExchanges.begin(),
		   parseSmallMolExchange(rMolUnit,
					 parsedPlex));

    // Parse the modification exchanges.
    xmlpp::Element* pModificationExchangesElt
      = domUtils::mustGetUniqueChild(pBndOmniGenElt,
				     eltName::modificationExchanges);
    xmlpp::Node::NodeList modificationExchangeNodes
      = pModificationExchangesElt->get_children(eltName::modificationExchange);
    std::vector<modificationExchange> modificationExchanges
      = std::vector<modificationExchange>(modificationExchangeNodes.size());
    std::transform(modificationExchangeNodes.begin(),
		   modificationExchangeNodes.end(),
		   modificationExchanges.begin(),
		   parseModificationExchange(rMolUnit,
					     parsedPlex));

    // Parse additional reactant.
    xmlpp::Element* pAdditionalReactantSpeciesElt
      = domUtils::getOptionalChild(pBndOmniGenElt,
				   eltName::additionalReactantSpecies);
    mzr::species* pAdditionalReactantSpecies = 0;
    if(pAdditionalReactantSpeciesElt)
      {
	std::string additionalReactantSpeciesName
	  = domUtils::mustGetAttrString
	  (pAdditionalReactantSpeciesElt,
	   eltName::additionalReactantSpecies_nameAttr);

	pAdditionalReactantSpecies
	  = rMzrUnit.mustFindSpecies(pAdditionalReactantSpeciesElt,
				     additionalReactantSpeciesName);
      }

    // Parse additional product.
    xmlpp::Element* pAdditionalProductSpeciesElt
      = domUtils::getOptionalChild(pBndOmniGenElt,
				   eltName::additionalProductSpecies);
    mzr::species* pAdditionalProductSpecies = 0;
    if(pAdditionalProductSpeciesElt)
      {
	std::string additionalProductSpeciesName
	  = domUtils::mustGetAttrString
	  (pAdditionalProductSpeciesElt,
	   eltName::additionalProductSpecies_nameAttr);

	pAdditionalProductSpecies
	  = rMzrUnit.mustFindSpecies(pAdditionalProductSpeciesElt,
				     additionalProductSpeciesName);
      }

    // Parse the reaction rate, whose units depend on whether an additional
    // reactant species is given or not.
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pBndOmniGenElt,
				     eltName::rate);
    double rate
      = domUtils::mustGetAttrPosDouble(pRateElt,
				       eltName::rate_valueAttr);

    // Construct the reaction rate exptrapolator.  The generated reactions
    // could be either unary or binary, depending on whether an additional
    // reactant is given, and there is a constructor for each of these cases.
    bndOmniMassExtrap* pExtrapolator = 0;
    if(pAdditionalReactantSpeciesElt)
      {
	// Construct the default species of the triggering omniplex.
	//
	// The need to construct a default member of a plexFamily arises
	// again.  Reinsitute it as a member variable of plexFamily?
	plx::plexFamily* pTriggeringFamily = pOmni->getFamily();
	std::vector<bnd::molParam> defaultParams
	  = pTriggeringFamily->makeDefaultMolParams();
	plx::plexSpecies* pDefaultTriggeringSpecies
	  = pTriggeringFamily->makeMember(defaultParams);

	// Make sure that the molecular weight of the additional reactant
	// species can be determined.
	mzr::massive* pAdditionalMassiveSpecies
	  = mzr::mustBeMassiveSpecies(pAdditionalReactantSpeciesElt,
				      pAdditionalReactantSpecies);

	pExtrapolator = new bndOmniMassExtrap(rate,
					      pDefaultTriggeringSpecies,
					      pAdditionalMassiveSpecies);
      }
    else
      {
	pExtrapolator = new bndOmniMassExtrap(rate);
      }

    // Construct the reaction family and register it for memory management.
    bndOmniFam* pFamily
      = new bndOmniFam(rMzrUnit,
		       rPlexUnit,
		       smallMolExchanges,
		       modificationExchanges,
		       pAdditionalReactantSpecies,
		       pAdditionalProductSpecies,
		       pExtrapolator);
    rMzrUnit.addReactionFamily(pFamily);

    // Connect the family's reaction generator to the triggering omniplex's
    // feature.
    pOmni->getSubPlexFeature()->addRxnGen(pFamily->getRxnGen());
  }
}
