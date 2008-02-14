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

#include "fnd/speciesNotMassiveXcpt.hh"
#include "cml/cptModMol.hh"
#include "cml/badModMolXcpt.hh"
#include "clx/parserPlex.hh"
#include "clx/parseOmniPlex.hh"
#include "cft/cftEltName.hh"
#include "cft/omniGen.hh"
#include "cft/omniFam.hh"
#include "cft/parseOmniGen.hh"

namespace cft
{
  class parseSmallMolExchange :
    public std::unary_function<xmlpp::Node*, smallMolExchange>
  {
    cml::cmlUnit& rMolUnit;
    const clx::parserPlex& rParsedPlex;
  public:
    parseSmallMolExchange(cml::cmlUnit& refMolUnit,
			  const clx::parserPlex& refParsedPlex) :
      rMolUnit(refMolUnit),
      rParsedPlex(refParsedPlex)
    {}

    smallMolExchange
    operator()(xmlpp::Node* pSmallMolExchangeNode) const
      throw(utl::xcpt)
    {
      xmlpp::Element* pSmallMolExchangeElt
	= utl::dom::mustBeElementPtr(pSmallMolExchangeNode);
      
      xmlpp::Element* pSmallMolInstanceRefElt
	= utl::dom::mustGetUniqueChild(pSmallMolExchangeElt,
				       eltName::smallMolInstanceRef);

      std::string smallMolInstanceName
	= utl::dom::mustGetAttrString(pSmallMolInstanceRefElt,
				      eltName::smallMolInstanceRef_nameAttr);
      cpx::molSpec exchangedMolSpec
	= rParsedPlex.mustGetMolNdxByName(pSmallMolInstanceRefElt,
					  smallMolInstanceName);

      xmlpp::Element* pSmallMolRefElt
	= utl::dom::mustGetUniqueChild(pSmallMolExchangeElt,
				       eltName::smallMolRef);

      std::string smallMolName
	= utl::dom::mustGetAttrString(pSmallMolRefElt,
				      eltName::smallMolRef_nameAttr);
      
      cml::cptSmallMol* pReplacementMol
	= rMolUnit.mustFindSmallMol(smallMolName,
				    pSmallMolInstanceRefElt);

      return smallMolExchange(exchangedMolSpec,
			      pReplacementMol);
    }
  };

  class parseModificationExchange :
    public std::unary_function<xmlpp::Node*, modificationExchange>
  {
    cml::cmlUnit& rMolUnit;
    const clx::parserPlex& rParsedPlex;
  public:
    parseModificationExchange(cml::cmlUnit& refMolUnit,
			      const clx::parserPlex& refParsedPlex) :
      rMolUnit(refMolUnit),
      rParsedPlex(refParsedPlex)
    {
    }

    modificationExchange
    operator()(xmlpp::Node* pModificationExchangeNode) const
      throw(utl::xcpt)
    {
      xmlpp::Element* pModificationExchangeElt
	= utl::dom::mustBeElementPtr(pModificationExchangeNode);

      xmlpp::Element* pModMolInstanceRefElt
	= utl::dom::mustGetUniqueChild(pModificationExchangeElt,
				       eltName::modMolInstanceRef);

      // Get the instance name of the mol instance in which the modification
      // exchange is to take place, and get the index of the mol in the
      // complex.
      std::string modMolInstanceName
	= utl::dom::mustGetAttrString(pModMolInstanceRefElt,
				      eltName::modMolInstanceRef_nameAttr);

      cpx::molSpec modMolSpec
	= rParsedPlex.mustGetMolNdxByName(pModMolInstanceRefElt,
					  modMolInstanceName);

      // Get the modification site name.
      xmlpp::Element* pModSiteRefElt
	= utl::dom::mustGetUniqueChild(pModMolInstanceRefElt,
				       eltName::modSiteRef);

      std::string modSiteName
	= utl::dom::mustGetAttrString(pModSiteRefElt,
				      eltName::modSiteRef_nameAttr);

      // Have to use the mol to translate modification site name to
      // modification index.
      const cml::cptModMol* pModMol
	= cml::mustBeModMol
	(rParsedPlex.mustGetMolByName(pModMolInstanceRefElt,
				      modMolInstanceName),
	 pModMolInstanceRefElt);

      int modSiteNdx
	= pModMol->mustGetModSiteNdx(modSiteName,
				     pModSiteRefElt);

      // Get the name of the modification that is to replace the modification
      // currently found at the modification site.
      xmlpp::Element* pInstalledModRefElt
	= utl::dom::mustGetUniqueChild(pModificationExchangeElt,
				       eltName::installedModRef);
      std::string modificationName
	= utl::dom::mustGetAttrString(pInstalledModRefElt,
				      eltName::installedModRef_nameAttr);

      // Look up the modification using the cmlUnit.
      const cpx::modification* pModification
	= rMolUnit.mustGetMod(modificationName,
			      pInstalledModRefElt);

      return modificationExchange(modMolSpec,
				  modSiteNdx,
				  pModification);
    }
  };

  void
  parseOmniGen::
  operator()(xmlpp::Node* pOmniGenNode) const
    throw(utl::xcpt)
  {
    xmlpp::Element* pOmniGenElt
      = utl::dom::mustBeElementPtr(pOmniGenNode);

    // Parse the enabling omniplex.
    xmlpp::Element* pEnablingOmniElt
      = utl::dom::mustGetUniqueChild(pOmniGenElt,
				     eltName::enablingOmniplex);
    clx::parserPlex parsedPlex;
    clx::cptOmniPlex* pOmni
      = clx::findOmni(pEnablingOmniElt,
		      rCmlUnit,
		      rClxUnit,
		      parsedPlex);

    // Parse the small-mol exchanges.
    xmlpp::Element* pSmallMolExchangesElt
      = utl::dom::mustGetUniqueChild(pOmniGenElt,
				     eltName::smallMolExchanges);
    xmlpp::Node::NodeList smallMolExchangeNodes
      = pSmallMolExchangesElt->get_children(eltName::smallMolExchange);
    std::vector<smallMolExchange> smallMolExchanges
      = std::vector<smallMolExchange>(smallMolExchangeNodes.size());
    std::transform(smallMolExchangeNodes.begin(),
		   smallMolExchangeNodes.end(),
		   smallMolExchanges.begin(),
		   parseSmallMolExchange(rCmlUnit,
					 parsedPlex));

    // Parse the modification exchanges.
    xmlpp::Element* pModificationExchangesElt
      = utl::dom::mustGetUniqueChild(pOmniGenElt,
				     eltName::modificationExchanges);
    xmlpp::Node::NodeList modificationExchangeNodes
      = pModificationExchangesElt->get_children(eltName::modificationExchange);
    std::vector<modificationExchange> modificationExchanges
      = std::vector<modificationExchange>(modificationExchangeNodes.size());
    std::transform(modificationExchangeNodes.begin(),
		   modificationExchangeNodes.end(),
		   modificationExchanges.begin(),
		   parseModificationExchange(rCmlUnit,
					     parsedPlex));

    // Parse additional reactant.
    xmlpp::Element* pAdditionalReactantSpeciesElt
      = utl::dom::getOptionalChild(pOmniGenElt,
				   eltName::additionalReactantSpecies);
    cpt::globalSpecies* pAdditionalReactantSpecies = 0;
    if(pAdditionalReactantSpeciesElt)
      {
	std::string additionalReactantSpeciesName
	  = utl::dom::mustGetAttrString
	  (pAdditionalReactantSpeciesElt,
	   eltName::additionalReactantSpecies_nameAttr);

	pAdditionalReactantSpecies
	  = rCptUnit.mustFindSpecies(additionalReactantSpeciesName,
				     pAdditionalReactantSpeciesElt);
      }

    // Parse additional product.
    xmlpp::Element* pAdditionalProductSpeciesElt
      = utl::dom::getOptionalChild(pOmniGenElt,
				   eltName::additionalProductSpecies);
    cpt::globalSpecies* pAdditionalProductSpecies = 0;
    if(pAdditionalProductSpeciesElt)
      {
	std::string additionalProductSpeciesName
	  = utl::dom::mustGetAttrString
	  (pAdditionalProductSpeciesElt,
	   eltName::additionalProductSpecies_nameAttr);

	pAdditionalProductSpecies
	  = rCptUnit.mustFindSpecies(additionalProductSpeciesName,
				     pAdditionalProductSpeciesElt);
      }

    // Parse the reaction rate, whose units depend on whether an additional
    // reactant species is given or not.
    xmlpp::Element* pRateElt
      = utl::dom::mustGetUniqueChild(pOmniGenElt,
				     eltName::rate);
    double rate
      = utl::dom::mustGetAttrPosDouble(pRateElt,
				       eltName::rate_valueAttr);

    // Construct the reaction rate exptrapolator.  The generated reactions
    // could be either unary or binary, depending on whether an additional
    // reactant is given, and there is a constructor for each of these cases.
    omniMassExtrap* pExtrapolator = 0;
    if(pAdditionalReactantSpeciesElt)
      {
	// Construct the default species of the triggering omniplex.
	//
	// The need to construct a default member of a plexFamily arises
	// again.  Reinsitute it as a member variable of plexFamily?
	clx::cptPlexFamily* pTriggeringFamily = pOmni->getFamily();
	std::vector<cpx::molParam> defaultParams
	  = pTriggeringFamily->makeDefaultMolParams();
	clx::cptPlexSpecies* pDefaultTriggeringSpecies
	  = pTriggeringFamily->makeMember(defaultParams);

	// Make sure that the molecular weight of the additional reactant
	// species can be determined.
	fnd::massive* pAdditionalMassiveSpecies
	  = fnd::mustBeMassiveSpecies(pAdditionalReactantSpecies,
				      pAdditionalReactantSpeciesElt);

	pExtrapolator = new omniMassExtrap(rate,
					   pDefaultTriggeringSpecies,
					   pAdditionalMassiveSpecies);
      }
    else
      {
	pExtrapolator = new omniMassExtrap(rate);
      }

    // Construct the reaction family and register it for memory management.
    omniFam* pFamily
      = new omniFam(rCptApp,
		    rCptUnit,
		    rClxUnit,
		    smallMolExchanges,
		    modificationExchanges,
		    pAdditionalReactantSpecies,
		    pAdditionalProductSpecies,
		    pExtrapolator);
    rCptUnit.addReactionFamily(pFamily);

    // Connect the family's reaction generator to the triggering omniplex's
    // feature.
    pOmni->getSubPlexFeature()->insert(pFamily->getRxnGen());
  }
}
