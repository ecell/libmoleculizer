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
#include "mol/modQuery.hh"
#include "ftr/ftrEltName.hh"
#include "ftr/uniMolGen.hh"
#include "ftr/uniMolFam.hh"
#include "ftr/parseUniMolGen.hh"

namespace ftr
{
  bnd::modMol*
  parseEnablingMol(xmlpp::Element* pEnablingMolElt,
		   bnd::molUnit& rMolUnit)
    throw(mzr::mzrXcpt)
  {
    std::string enablingMolName
      = domUtils::mustGetAttrString(pEnablingMolElt,
				    eltName::enablingMol_nameAttr);

    return rMolUnit.mustFindModMol(pEnablingMolElt,
				   enablingMolName);
  }

  class parseEnablingModSiteRef :
    public std::unary_function<xmlpp::Node*, void>
  {
    bnd::modMol* pEnablingMol;
    bnd::andMolQueries* pQuery;
    bnd::molUnit& rMolUnit;
  public:
    parseEnablingModSiteRef(bnd::modMol* pEnablingModMol,
			    bnd::andMolQueries* pEnableQuery,
			    bnd::molUnit& refMolUnit) :
      pEnablingMol(pEnablingModMol),
      pQuery(pEnableQuery),
      rMolUnit(refMolUnit)
    {}

    void
    operator()(xmlpp::Node* pModSiteRefNode) const
      throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pModSiteRefElt
	= domUtils::mustBeElementPtr(pModSiteRefNode);

      // Parse the modification site name and convert it into an index.
      std::string modSiteName
	= domUtils::mustGetAttrString(pModSiteRefElt,
				      eltName::modSiteRef_nameAttr);
      int modSiteNdx
	= pEnablingMol->mustGetModSiteNdx(pModSiteRefElt,
					  modSiteName);

      // Parse the modifcation name, and get the modification.
      xmlpp::Element* pModRefElt
	= domUtils::mustGetUniqueChild(pModSiteRefElt,
				       eltName::modRef);
      std::string modName
	= domUtils::mustGetAttrString(pModRefElt,
				      eltName::modRef_nameAttr);
      const bnd::modification* pModToSee
	= rMolUnit.mustGetMod(pModRefElt,
			      modName);

      // Add the test that the modification site have the given modification
      // to the conjunction of mol state queries.
      pQuery->addQuery(new bnd::modMixinMolQuery(rMolUnit,
						 modSiteNdx,
						 pModToSee));
    }
  };
  
  bnd::andMolQueries*
  parseEnablingModifications(xmlpp::Element* pEnablingModsElt,
			     bnd::modMol* pEnablingMol,
			     bnd::molUnit& rMolUnit)
    throw(mzr::mzrXcpt)
  {
    bnd::andMolQueries* pQuery
      = new bnd::andMolQueries(rMolUnit);

    xmlpp::Node::NodeList modSiteRefNodes
      = pEnablingModsElt->get_children(eltName::modSiteRef);

    std::for_each(modSiteRefNodes.begin(),
		  modSiteRefNodes.end(),
		  parseEnablingModSiteRef(pEnablingMol,
					  pQuery,
					  rMolUnit));
    return pQuery;
  }

  class parseModificationExchange :
    public std::unary_function<xmlpp::Node*, molModExchange>
  {
    const bnd::molUnit& rMolUnit;
    const bnd::modMol* pEnabling;
  public:
    parseModificationExchange(const bnd::molUnit& refMolUnit,
			      const bnd::modMol* pEnablingModMol) :
      rMolUnit(refMolUnit),
      pEnabling(pEnablingModMol)
    {}

    molModExchange
    operator()(xmlpp::Node* pModificationExchangeNode) const
      throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pModificationExchangeElt
	= domUtils::mustBeElementPtr(pModificationExchangeNode);

      // Parse the name of the modification site to be changed, and
      // convert it to site index.
      xmlpp::Element* pModSiteRefElt
	= domUtils::mustGetUniqueChild(pModificationExchangeElt,
				       eltName::modSiteRef);
      std::string modSiteName
	= domUtils::mustGetAttrString(pModSiteRefElt,
				      eltName::modSiteRef_nameAttr);
      int modSiteNdx
	= pEnabling->mustGetModSiteNdx(pModSiteRefElt,
				       modSiteName);

      // Parse the name of the modification to be substituted in, and
      // look up the corresponding modification.
      xmlpp::Element* pInstalledModRefElt
	= domUtils::mustGetUniqueChild(pModificationExchangeElt,
				       eltName::installedModRef);
      std::string modName
	= domUtils::mustGetAttrString(pInstalledModRefElt,
				      eltName::installedModRef_nameAttr);
      const bnd::modification* pMod
	= rMolUnit.mustGetMod(pInstalledModRefElt,
			      modName);

      return molModExchange(modSiteNdx,
			    pMod);
    }
  };

  std::vector<molModExchange>
  parseModificationExchanges(xmlpp::Element* pModificationExchangesElt,
			     const bnd::modMol* pEnablingModMol,
			     const bnd::molUnit& rMolUnit)
    throw(mzr::mzrXcpt)
  {
    
    xmlpp::Node::NodeList modificationExchangeNodes
      = pModificationExchangesElt->get_children(eltName::modificationExchange);

    std::vector<molModExchange> modificationExchanges
      = std::vector<molModExchange>(modificationExchangeNodes.size());

    std::transform(modificationExchangeNodes.begin(),
		   modificationExchangeNodes.end(),
		   modificationExchanges.begin(),
		   parseModificationExchange(rMolUnit,
					     pEnablingModMol));

    return modificationExchanges;
  }

  mzr::species*
  parseAdditionalReactantSpecies(xmlpp::Element* pAdditionalReactantSpeciesElt,
				 const mzr::mzrUnit& rMzrUnit)
    throw(mzr::mzrXcpt)
  {
    std::string speciesName
      = domUtils::mustGetAttrString
      (pAdditionalReactantSpeciesElt,
       eltName::additionalReactantSpecies_nameAttr);

    return rMzrUnit.mustFindSpecies(pAdditionalReactantSpeciesElt,
				    speciesName);
  }

  mzr::species*
  parseAdditionalProductSpecies(xmlpp::Element* pAdditionalProductSpeciesElt,
				const mzr::mzrUnit& rMzrUnit)
    throw(mzr::mzrXcpt)
  {
    std::string speciesName
      = domUtils::mustGetAttrString
      (pAdditionalProductSpeciesElt,
       eltName::additionalProductSpecies_nameAttr);

    return rMzrUnit.mustFindSpecies(pAdditionalProductSpeciesElt,
				    speciesName);
  }

  void
  parseUniMolGen::
  operator()(xmlpp::Node* pUniMolGenNode) const
    throw(mzr::mzrXcpt)
  {
    xmlpp::Element* pUniMolGenElt
      = domUtils::mustBeElementPtr(pUniMolGenNode);

    // Parse the enabling mol.
    xmlpp::Element* pEnablingMolElt
      = domUtils::mustGetUniqueChild(pUniMolGenElt,
				     eltName::enablingMol);
    bnd::modMol* pEnablingModMol
      = parseEnablingMol(pEnablingMolElt,
			 rMolUnit);

    // Parse the enabling modifications.
    xmlpp::Element* pEnablingModsElt
      = domUtils::mustGetUniqueChild(pUniMolGenElt,
				     eltName::enablingModifications);
    bnd::andMolQueries* pQuery
      = parseEnablingModifications(pEnablingModsElt,
				   pEnablingModMol,
				   rMolUnit);

    // Parse the modification exchanges.
    xmlpp::Element* pModificationExchangesElt
      = domUtils::mustGetUniqueChild(pUniMolGenElt,
				     eltName::modificationExchanges);

    std::vector<molModExchange> modificationExchanges
      = parseModificationExchanges(pModificationExchangesElt,
				   pEnablingModMol,
				   rMolUnit);

    // Parse the reation rate.
    //
    // (This is out of order, so that we can use the rate in the construction
    // of the rate extrapolator.  Which rate extrapolator is used depends
    // on whether there is an additional reactant.)
    xmlpp::Element* pRateElt
      = domUtils::mustGetUniqueChild(pUniMolGenElt,
				     eltName::rate);
    double rate
      = domUtils::mustGetAttrDouble(pRateElt,
				    eltName::rate_valueAttr);

    // Parse the additional reactant species, if any.
    // 
    // If there is one, use the binary reaction rate extrapolator,
    // which requires that the additional reactant be massive.
    // If there isn't one, use the unary reaction rate extrapolator,
    // which doesn't really do anything.
    xmlpp::Element* pAdditionalReactantSpeciesElt
      = domUtils::getOptionalChild(pUniMolGenElt,
				   eltName::additionalReactantSpecies);
    
    mzr::species* pAdditionalReactantSpecies = 0;
    uniMolExtrapolator* pExtrapolator = 0;
    if(pAdditionalReactantSpeciesElt)
      {
	pAdditionalReactantSpecies
	  = parseAdditionalReactantSpecies(pAdditionalReactantSpeciesElt,
					   rMzrUnit);

	mzr::massive* pMassive
	  = mzr::mustBeMassiveSpecies(pAdditionalReactantSpeciesElt,
				      pAdditionalReactantSpecies);

	pExtrapolator = new uniMolMassExtrap(rate,
					     pEnablingModMol,
					     pMassive);
      }
    else
      {
	pExtrapolator = new uniMolMassExtrap(rate);
      }

    // Parse the additional product species, if any.
    xmlpp::Element* pAdditionalProductSpeciesElt
      = domUtils::getOptionalChild(pUniMolGenElt,
				   eltName::additionalProductSpecies);

    mzr::species* pAdditionalProductSpecies = 0;
    if(pAdditionalProductSpeciesElt)
      {
	pAdditionalProductSpecies
	  = parseAdditionalProductSpecies(pAdditionalProductSpeciesElt,
					  rMzrUnit);
      }

    // Construct the reaction family and register it for memory management.
    uniMolFam* pFamily
      = new uniMolFam(rMzrUnit,
		      rPlexUnit,
		      pEnablingModMol,
		      pQuery,
		      modificationExchanges,
		      pAdditionalReactantSpecies,
		      pAdditionalProductSpecies,
		      pExtrapolator);
    rMzrUnit.addReactionFamily(pFamily);

    // Connect the family's reaction generator to the mol's feature.
    // Note that mol inherits from feature.
    pEnablingModMol->addRxnGen(pFamily->getRxnGen());
  }
}
