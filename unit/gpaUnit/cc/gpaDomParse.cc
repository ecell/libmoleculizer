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
#include "mzr/mzrEltName.hh"
#include "mol/molDomParse.hh"
#include "mol/molUnit.hh"
#include "mol/molEltName.hh"
#include "plex/plexDomParse.hh"
#include "plex/plexFamily.hh"
#include "plex/plexEltName.hh"
#include "stoch/stochUnit.hh"
#include "gpa/gpaExchangeFam.hh"
#include "gpa/gpaRevertFam.hh"
#include "gpa/gpaUnit.hh"
#include "gpa/gpaEltName.hh"

namespace gpa
{
  // The complex identified as the "enabling complex" for the GpaExchange
  // reaction needs to be identified by the plex unit as an omniplex
  // (for recognition of all species of complex where the GpaExchange
  // reaction can take place) and as a particular species, in order to
  // get its mass and calculate the weight-invariant reaction rate constant.
  //
  // This function produces the particular species of complex, assuming that
  // the plex unit has previously scanned and found the complex as an omniplex.
  class addGpaExchangeRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    stoch::stochUnit& rStochUnit;
    
  public:
    addGpaExchangeRxnFam(mzr::mzrUnit& refMzrUnit,
			 bnd::molUnit& refMolUnit,
			 plx::plexUnit& refPlexUnit,
			 stoch::stochUnit& refStochUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      rStochUnit(refStochUnit)
    {}
    
    void
    operator()(xmlpp::Node* pGpaExchangeGenNode) const throw(std::exception)
    {
      xmlpp::Element* pGpaExchangeGenElt
	= domUtils::mustBeElementPtr(pGpaExchangeGenNode);

      // Extract the enabling plex's element.
      //
      // Note that the plexUnit will have to process this plex
      // in its initial pass.  This means, in general, that the plex
      // unit will have code in it that depends on all the other units
      // that have to do this.
      //
      // Maybe I can arrange a more general way to find all of these
      // for the plexUnit, perhaps make sure that they are all "plex"
      // elements?  Need something analogous for omniplexes....
      xmlpp::Element* pEnablingPlexElt
	= domUtils::mustGetUniqueChild(pGpaExchangeGenElt,
				       plx::eltName::plex);

      // Parse the enabling complex.  The plex unit will also have
      // to be able to find this complex....
      plx::parserPlex enablingPlex;
      plx::plexFamily* pEnablingPlexFamily
	= recognizePlexElt(pEnablingPlexElt,
			   enablingPlex,
			   rMolUnit,
			   rPlexUnit);

      // Determine the index of the gpa mol in the enabling complex.
      xmlpp::Element* pTargetModMolInstanceRefElt
	= domUtils::mustGetUniqueChild(pGpaExchangeGenElt,
				       eltName::targetModMolInstanceRef);
      std::string modMolInstanceName
	= domUtils::mustGetAttrString
	(pTargetModMolInstanceRefElt,
	 eltName::targetModMolInstanceRef_nameAttr);

      // Get the gpa mol as a basic mol, as well as its index in the plex.
      //
      // This exception is in plexDomParse.hh.
      int gpaMolSpec = enablingPlex.getMolNdxByName(modMolInstanceName);
      if(gpaMolSpec < 0)
	throw plx::unknownMolInstanceXcpt(pTargetModMolInstanceRefElt,
					  modMolInstanceName);
      bnd::mol* pGpaMol = enablingPlex.mols[gpaMolSpec];

      // Make sure the gpa mol is a mod-mol.
      //
      // This exception is in molDomParse.hh
      bnd::modMol* pGpaModMol = dynamic_cast<bnd::modMol*>(pGpaMol);
      if(0 == pGpaModMol) throw bnd::badModMolCastXcpt(pTargetModMolInstanceRefElt,
						  pGpaMol);

      // Get the name of the gtp binding site on the gpa mol.
      xmlpp::Element* pModSiteRefElt
	= domUtils::mustGetUniqueChild(pTargetModMolInstanceRefElt,
				       bnd::eltName::modSiteRef);
      std::string modSiteName
	= domUtils::mustGetAttrString(pModSiteRefElt,
				      bnd::eltName::modSiteRef_nameAttr);

      // Convert the modification site name into modification site index.
      //
      // This exception is in molDomParse.hh.
      int gtpModSiteNdx = pGpaModMol->getModSiteNdx(modSiteName);
      if(gtpModSiteNdx < 0) throw bnd::unknownModSiteXcpt(pModSiteRefElt,
						     modSiteName,
						     pGpaModMol);

      // Get the "gdp-bound" modification.
      //
      // molUnit::mustGetMod is in molUnit.hh.
      xmlpp::Element* pGdpBoundModRefElt
	= domUtils::mustGetUniqueChild(pGpaExchangeGenElt,
				       eltName::gdpBoundModRef);
      std::string gdpBoundModName
	= domUtils::mustGetAttrString(pGdpBoundModRefElt,
				      eltName::gdpBoundModRef_nameAttr);
      const bnd::modification* pGdpBoundMod
	= rMolUnit.mustGetMod(gdpBoundModName);

      // Get the "gtp-bound" modification.
      xmlpp::Element* pGtpBoundModRefElt
	= domUtils::mustGetUniqueChild(pGpaExchangeGenElt,
				       eltName::gtpBoundModRef);
      std::string gtpBoundModName
	= domUtils::mustGetAttrString(pGtpBoundModRefElt,
				      eltName::gtpBoundModRef_nameAttr);
      const bnd::modification* pGtpBoundMod
	= rMolUnit.mustGetMod(gtpBoundModName);

      // Get the GDP species.
      //
      // unknownSpeciesXcpt is in mzrUnitParse.hh.
      xmlpp::Element* pGdpSpeciesElt
	= domUtils::mustGetUniqueChild(pGpaExchangeGenElt,
				       eltName::gdpSpeciesRef);
      std::string gdpSpeciesName
	= domUtils::mustGetAttrString(pGdpSpeciesElt,
				      eltName::gdpSpeciesRef_nameAttr);

      stoch::stochSpecies* pGDP
	= rStochUnit.mustGetStochSpecies(pGdpSpeciesElt,
					 gdpSpeciesName);

      // Get the GTP species.
      xmlpp::Element* pGtpSpeciesElt
	= domUtils::mustGetUniqueChild(pGpaExchangeGenElt,
				       eltName::gtpSpeciesRef);
      std::string gtpSpeciesName
	= domUtils::mustGetAttrString(pGtpSpeciesElt,
				      eltName::gtpSpeciesRef_nameAttr);

      stoch::stochSpecies* pGTP
	= rStochUnit.mustGetStochSpecies(pGtpSpeciesElt,
					 gtpSpeciesName);

      // Construct the plex species, borrowing a routine from plexDomParse.cc.
      //
      // In fact, this bit of code is modified form of "processPlexSpecies"
      // from plexDomParse.cc.
      std::vector<bnd::molParam> molParams
	= pEnablingPlexFamily->makeDefaultMolParams();
      xmlpp::Element* pInstanceStatesElt
	= domUtils::mustGetUniqueChild(pGpaExchangeGenElt,
				       plx::eltName::instanceStates);
      xmlpp::Node::NodeList modMolInstanceRefs
	= pInstanceStatesElt->get_children(plx::eltName::modMolInstanceRef);
      std::for_each(modMolInstanceRefs.begin(),
		    modMolInstanceRefs.end(),
		    plx::replaceModMolDefaultState(molParams,
						   pEnablingPlexFamily,
						   enablingPlex,
						   rMolUnit));

      // This is a misnomer; the subcomplex is enabling, the species is just
      // for weight calculation.
      plx::plexSpecies* pEnablingSpecies
	= pEnablingPlexFamily->getMember(molParams);

      // Get the rate, and combine it with the masses of GTP and the
      // "enabling species" to get the weight-invariant form
      // of the rate constant.
      xmlpp::Element* pRateElt
	= domUtils::mustGetUniqueChild(pGpaExchangeGenElt,
				       mzr::eltName::rate);
      double rate
	= domUtils::mustGetAttrDouble(pRateElt,
				      mzr::eltName::rate_valueAttr);

      // Compute the weight-invariant form of the rate,
      // using the molecular weight of the complex and the molecular weight
      // of GTP.
      //
      // But for this, we have to realize both the complex and GTP as
      // massive species, which just adds another complication to the
      // complication with the plexUnit's having to find this
      // enabling complex.
      //
      // We have to construct a particular
//       double rateConstant = mzr::bindingInvariant(rate,
// 						  pGTP->getWeight(),
// 						  pEnablingSpecies->getWeight());

      // Construct reaction rate extrapolator.  What kind of extrapolator
      // is constructed here will soon be user-selectable.  Extrapolators
      // are deleted by the reaction generators they are attached to.
      gpaExchangeMassExtrap* pExtrap
	= new gpaExchangeMassExtrap(rate,
				    pEnablingSpecies,
				    pGTP);

      // Construct the gpa exchange reaction family.
      gpaExchangeFam* pReactionFamily
	= new gpaExchangeFam(pGpaModMol,
			     gpaMolSpec,
			     gtpModSiteNdx,
			     pGdpBoundMod,
			     pGtpBoundMod,
			     pGDP,
			     pGTP,
			     rMzrUnit,
			     pExtrap);

      rMzrUnit.addReactionFamily(pReactionFamily);

      // Connect the family's reaction generator to the omniplex
      // feature of the enabling complex.
      plx::omniPlexFeature* pEnablingOmniPlexFeature
	= pEnablingPlexFamily->getSubPlexFeature();
      pEnablingOmniPlexFeature->addRxnGen(pReactionFamily->getRxnGen());
    }
  };

  class addGpaRevertRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    
  public:
    addGpaRevertRxnFam(mzr::mzrUnit& refMzrUnit,
		       bnd::molUnit& refMolUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit)
    {}
    
    void
    operator()(xmlpp::Node* pGpaRevertGenNode) const throw(std::exception)
    {
      xmlpp::Element* pGpaRevertGenElt
	= domUtils::mustBeElementPtr(pGpaRevertGenNode);

      // Get the gpa modMol.
      xmlpp::Element* pModMolRefElt
	= domUtils::mustGetUniqueChild(pGpaRevertGenElt,
				       eltName::modMolRef);
      std::string gpaMolName
	= domUtils::mustGetAttrString(pModMolRefElt,
				      eltName::modMolRef_nameAttr);

      bnd::mol* pGpaMol = rMolUnit.findMol(gpaMolName);
      if(0 == pGpaMol) throw bnd::unknownMolXcpt(pModMolRefElt,
					    gpaMolName);

      // This exception is in molDomParse.hh.
      bnd::modMol* pGpaModMol = dynamic_cast<bnd::modMol*>(pGpaMol);
      if(0 == pGpaModMol) throw bnd::badModMolCastXcpt(pModMolRefElt,
						  pGpaMol);

      // Get the modification site at which GTP/GDP bind.
      xmlpp::Element* pModSiteRefElt
	= domUtils::mustGetUniqueChild(pModMolRefElt,
				       bnd::eltName::modSiteRef);

      std::string modSiteName
	= domUtils::mustGetAttrString(pModSiteRefElt,
				      bnd::eltName::modSiteRef_nameAttr);
      int gtpModSiteNdx
	= pGpaModMol->getModSiteNdx(modSiteName);
      if(gtpModSiteNdx < 0) throw bnd::unknownModSiteXcpt(pModSiteRefElt,
						     modSiteName,
						     pGpaModMol);

      // Get pointer to the gtp-bound modification.
      xmlpp::Element* pGtpBoundModRefElt
	= domUtils::mustGetUniqueChild(pGpaRevertGenElt,
				       eltName::gtpBoundModRef);
      std::string gtpBoundModName
	= domUtils::mustGetAttrString(pGtpBoundModRefElt,
				      eltName::gtpBoundModRef_nameAttr);
      const bnd::modification* pGtpBoundMod
	= rMolUnit.getMod(gtpBoundModName);
      if(0 == pGtpBoundMod) throw bnd::unknownModXcpt(pGtpBoundModRefElt,
						      gtpBoundModName);

      // Get pointer to the gdp-bound modification.
      xmlpp::Element* pGdpBoundModRefElt
	= domUtils::mustGetUniqueChild(pGpaRevertGenElt,
				       eltName::gdpBoundModRef);
      std::string gdpBoundModName
	= domUtils::mustGetAttrString(pGdpBoundModRefElt,
				      eltName::gdpBoundModRef_nameAttr);
      const bnd::modification* pGdpBoundMod
	= rMolUnit.getMod(gdpBoundModName);
      if(0 == pGdpBoundMod) throw bnd::unknownModXcpt(pGdpBoundModRefElt,
						 gdpBoundModName);

      // Get pointer to the phosphate species.
      xmlpp::Element* pPhosphateSpeciesRefElt
	= domUtils::mustGetUniqueChild(pGpaRevertGenElt,
				       eltName::phosphateSpeciesRef);

      std::string phosphateSpeciesName
	= domUtils::mustGetAttrString
	(pPhosphateSpeciesRefElt,
	 eltName::phosphateSpeciesRef_nameAttr);
    
      mzr::species* pPhosphateSpecies
	= rMzrUnit.findSpecies(phosphateSpeciesName);
      if(0 == pPhosphateSpecies)
	throw mzr::unknownSpeciesXcpt(pPhosphateSpeciesRefElt,
				 phosphateSpeciesName);

      // Get the reaction rate.
      xmlpp::Element* pRateElt
	= domUtils::mustGetUniqueChild(pGpaRevertGenElt,
				       mzr::eltName::rate);
      double rate
	= domUtils::mustGetAttrDouble(pRateElt,
				      mzr::eltName::rate_valueAttr);

      // Construct the reaction rate extrapolator.  Which extrapolator
      // is constructed will eventually be user selectable.  The extrapolator
      // is memory-managed by the reaction generator that it is attached to.
      gpaRevertNoExtrap* pExtrap
	= new gpaRevertNoExtrap(rate);

      // Construct the gpa revert reaction family.
      gpaRevertFam* pRevertFam
	= new gpaRevertFam(pGpaModMol,
			   gtpModSiteNdx,
			   pGdpBoundMod,
			   pGtpBoundMod,
			   pPhosphateSpecies,
			   rMzrUnit,
			   pExtrap);

      rMzrUnit.addReactionFamily(pRevertFam);

      // Add the reaction generator to the mol's list, where the mol is playing
      // its base class role as feature<plexSpecies, plexMolSpec>.
      pGpaModMol->addRxnGen(pRevertFam->getRxnGen());
    }
  };

  void
  gpaUnit::parseDomInput(xmlpp::Element* pRootElement,
			 xmlpp::Element* pModelElement,
			 xmlpp::Element* pStreamsElement,
			 xmlpp::Element* pEventsElement) throw(std::exception)
  {
    // Get the header node for all reaction generators.
    xmlpp::Element* pReactionGensElt
      = domUtils::mustGetUniqueChild(pModelElement,
				     mzr::eltName::reactionGens);

    // Get the gpa-exchange-gen elements.
    xmlpp::Node::NodeList gpaExchangeGenNodes
      = pReactionGensElt->get_children(eltName::gpaExchangeGen);

    // Add a gpa exchange reaction family for each element.
    std::for_each(gpaExchangeGenNodes.begin(),
		  gpaExchangeGenNodes.end(),
		  addGpaExchangeRxnFam(rMzrUnit,
				       rMolUnit,
				       rPlexUnit,
				       rStochUnit));

    // Get the gpa-revert-gen elements.
    xmlpp::Node::NodeList gpaRevertGenNodes
      = pReactionGensElt->get_children(eltName::gpaRevertGen);

    // Add a gpa revert reaction family for each element.
    std::for_each(gpaRevertGenNodes.begin(),
		  gpaRevertGenNodes.end(),
		  addGpaRevertRxnFam(rMzrUnit,
				     rMolUnit));
  }
}
