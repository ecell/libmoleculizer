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

#include <iostream> // For debugging.
#include <vector>
#include "mzr/pchem.hh"
#include "mzr/moleculizer.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/mzrEltName.hh"
#include "mol/molDomParse.hh"
#include "plex/plexUnit.hh"
#include "plex/plexEltName.hh"
#include "dimer/dimerUnit.hh"
#include "dimer/dimerizeFam.hh"
#include "dimer/decompFam.hh"
#include "dimer/dimerEltName.hh"

namespace dimer
{
  // I may need an anonymous namespace for this.  Should be doing that
  // in many places.
  typedef std::pair<bnd::mol*, int> bindingPartner;

  class parseBindingPartner :
    public std::unary_function<xmlpp::Node*, bindingPartner>
  {
    bnd::molUnit& rMolUnit;
    
  public:
    parseBindingPartner(bnd::molUnit& refMolUnit) :
      rMolUnit(refMolUnit)
    {}
    
    bindingPartner
    operator()(xmlpp::Node* pMolRefNode)
    {
      // Cast the mol-ref node to element; probably unnecessarily dynamically.
      xmlpp::Element* pMolRefElt
	= domUtils::mustBeElementPtr(pMolRefNode);

      // Get the mol name.
      std::string molName
	= domUtils::mustGetAttrString(pMolRefElt,
				      plx::eltName::molRef_nameAttr);

      // Why is the BAD BAD BAD static catalog of mols in plexUnit?
      bnd::mol* pMol = rMolUnit.findMol(molName);
      if(0 == pMol) throw bnd::unknownMolXcpt(pMolRefElt,
					      molName);

      // Get the binding site name.
      xmlpp::Element* pSiteRefElt
	= domUtils::mustGetUniqueChild(pMolRefElt,
				       eltName::siteRef);
      std::string siteName
	= domUtils::mustGetAttrString(pSiteRefElt,
				      eltName::siteRef_nameAttr);

      // Ask the mol to look up the index of the binding site
      // using the binding site name.
      int siteNdx = -1;
      if(! pMol->findSite(siteName,
			  siteNdx))
	throw bnd::unknownSiteXcpt(pSiteRefElt,
				   siteName);

      return std::make_pair(pMol,
			    siteNdx);
    }
  };

  class installDefaultKinetics : public
  std::unary_function<const std::pair<const std::string, bnd::siteShape>&, void>
  {
    const bnd::siteShape* pLeftShape;
    double onR;
    double offR;

    // Debugging matter.
    const bnd::mol* pLMol;
    const bnd::bindingSite& rLSite;
    const bnd::mol* pRMol;
    const bnd::bindingSite& rRSite;

    dimerizeExtrapolator* pDimerizeExtrap;
    decomposeExtrapolator* pDecomposeExtrap;
  
  public:
    installDefaultKinetics(const bnd::siteShape* pLeftSiteShape,
			   double onRate,
			   double offRate,
			   const bnd::mol* pLeftMol,
			   const bnd::bindingSite& rLeftSite,
			   const bnd::mol* pRightMol,
			   const bnd::bindingSite& rRightSite,
			   dimerizeExtrapolator* pDimerizeExtrapolator,
			   decomposeExtrapolator* pDecomposeExtrapolator) :
      pLeftShape(pLeftSiteShape),
      onR(onRate),
      offR(offRate),
      pLMol(pLeftMol),
      rLSite(rLeftSite),
      pRMol(pRightMol),
      rRSite(rRightSite),
      pDimerizeExtrap(pDimerizeExtrapolator),
      pDecomposeExtrap(pDecomposeExtrapolator)
    {}

    // The 'const' on the string key here DOES make a difference.
    // Without it, the compiler generates a local pair, apparently
    // because it sees a type difference.  It is enabled to do
    // something about it by the const'ness of the pair argument,
    // which legitimizes the construction of a temporary pair
    // of the right type, which again it seems to be able to do.
    //
    // Never mind that the pair is const, so that the string key is implicitly
    // const: const pair<int, int> is apparently not the same as const
    // pair<const int, const int>.
    void
    operator()(const std::pair<const std::string, bnd::siteShape>&
	       rRightSiteShapeEntry)
      const throw(mzr::mzrXcpt)
    {
      const bnd::siteShape* pRightShape = &(rRightSiteShapeEntry.second);

      pDimerizeExtrap->setRate(pLeftShape,
			       pRightShape,
			       onR);

      pDecomposeExtrap->setRate(pLeftShape,
				pRightShape,
				offR);
    }
  };

  class installDefaultKineticsToLeft : public
  std::unary_function<const std::pair<const std::string, bnd::siteShape>&, void>
  {
    const std::map<std::string, bnd::siteShape>& rRightShapeMap;
    double onR;
    double offR;

    // Debugging matter.
    const bnd::mol* pLMol;
    const bnd::bindingSite& rLSite;
    const bnd::mol* pRMol;
    const bnd::bindingSite& rRSite;

    dimerizeExtrapolator* pDimerizeExtrap;
    decomposeExtrapolator* pDecomposeExtrap;
  
  public:
    // Note that most of the arguments here are for debugging.
    installDefaultKineticsToLeft
    (const std::map<std::string, bnd::siteShape>& rRightSiteShapeMap,
     double onRate,
     double offRate,
     const bnd::mol* pLeftMol,
     const bnd::bindingSite& rLeftSite,
     const bnd::mol* pRightMol,
     const bnd::bindingSite& rRightSite,
     dimerizeExtrapolator* pDimerizeExtrapolator,
     decomposeExtrapolator* pDecomposeExtrapolator) :
      rRightShapeMap(rRightSiteShapeMap),
      onR(onRate),
      offR(offRate),
      pLMol(pLeftMol),
      rLSite(rLeftSite),
      pRMol(pRightMol),
      rRSite(rRightSite),
      pDimerizeExtrap(pDimerizeExtrapolator),
      pDecomposeExtrap(pDecomposeExtrapolator)
    {}

    void
    operator()(const std::pair<const std::string, bnd::siteShape>&
	       rLeftSiteShapeEntry)
      const throw(mzr::mzrXcpt)
    {
      std::for_each(rRightShapeMap.begin(),
		    rRightShapeMap.end(),
		    installDefaultKinetics(&(rLeftSiteShapeEntry.second),
					   onR,
					   offR,
					   pLMol,
					   rLSite,
					   pRMol,
					   rRSite,
					   pDimerizeExtrap,
					   pDecomposeExtrap));
    }
  };

  class getNameAttr :
    public std::unary_function<xmlpp::Node*, std::string>
  {
  public:
    std::string
    operator()(xmlpp::Node* pSiteShapeRefNode) const throw(mzr::mzrXcpt)
    {
      // Cast the site-shape-ref node* to element*; probably unnecessarily
      // dynamically.
      xmlpp::Element* pSiteShapeRefElt
	= domUtils::mustBeElementPtr(pSiteShapeRefNode);

      // Get the name of the site shape.
      return domUtils::mustGetAttrString(pSiteShapeRefElt,
					 bnd::eltName::siteShapeRef_nameAttr);
    }
  };

  class insertAlloRates :
    public std::unary_function<xmlpp::Node*, void>
  {
    const bnd::bindingSite& rLeftSite;
    double leftWeight;
  
    const bnd::bindingSite& rRightSite;
    double rightWeight;

    dimerizeExtrapolator* pDimerizeExtrap;
    decomposeExtrapolator* pDecomposeExtrap;
  
  public:
    insertAlloRates(const bnd::bindingSite& rLeftBindingSite,
		    double leftMolWeight,
		    const bnd::bindingSite& rRightBindingSite,
		    double rightMolWeight,
		    dimerizeExtrapolator* pDimerizeExtrapolator,
		    decomposeExtrapolator* pDecomposeExtrapolator) :
      rLeftSite(rLeftBindingSite),
      leftWeight(leftMolWeight),
      rRightSite(rRightBindingSite),
      rightWeight(rightMolWeight),
      pDimerizeExtrap(pDimerizeExtrapolator),
      pDecomposeExtrap(pDecomposeExtrapolator)
    {}

    void
    operator()(xmlpp::Node* pAlloRatesNode) const throw(mzr::mzrXcpt)
    {
      // Cast the allo-rates node* to element*; probably unnecessarily
      // dynamically.
      xmlpp::Element* pAlloRatesElt
	= domUtils::mustBeElementPtr(pAlloRatesNode);

      // Parse the site-shape-ref child elements, of which there must be 2,
      // into a more easily accessible form.
      xmlpp::Node::NodeList siteShapeRefNodes
	= pAlloRatesElt->get_children(bnd::eltName::siteShapeRef);
      if(siteShapeRefNodes.size() != 2)
	throw domUtils::badChildCountXcpt(pAlloRatesNode,
					  bnd::eltName::siteShapeRef,
					  2,
					  (int) siteShapeRefNodes.size());
      std::vector<std::string> siteShapeNames(2);
      std::transform(siteShapeRefNodes.begin(),
		     siteShapeRefNodes.end(),
		     siteShapeNames.begin(),
		     getNameAttr());

      // Get the corresponding site shape pointers.
      bnd::siteParam leftSiteParam
	= rLeftSite.getParam(siteShapeNames[0]);
      bnd::siteParam rightSiteParam
	= rRightSite.getParam(siteShapeNames[1]);

      // Get the on-rate.
      xmlpp::Element* pOnRateElt
	= domUtils::mustGetUniqueChild(pAlloRatesElt,
				       eltName::onRate);
      double onRate
	= domUtils::mustGetAttrDouble(pOnRateElt,
				      eltName::onRate_valueAttr);

      // Get the off-rate.
      xmlpp::Element* pOffRateElt
	= domUtils::mustGetUniqueChild(pAlloRatesElt,
				       eltName::offRate);
      double offRate
	= domUtils::mustGetAttrDouble(pOffRateElt,
				      eltName::offRate_valueAttr);

      // Install the allosteric rates, overwriting the default entries
      // which are already there.
      pDimerizeExtrap->setRate(leftSiteParam,
			       rightSiteParam,
			       onRate);

      pDecomposeExtrap->setRate(leftSiteParam,
				rightSiteParam,
				offRate);
    }
  };

  class addDimerizeDecompose :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    dimerUnit& rDimerUnit;
    plx::plexUnit& rPlexUnit;
    
  public:
    addDimerizeDecompose(mzr::mzrUnit& refMzrUnit,
			 bnd::molUnit& refMolUnit,
			 dimerUnit& refDimerUnit,
			 plx::plexUnit& refPlexUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rDimerUnit(refDimerUnit),
      rPlexUnit(refPlexUnit)
    {}
    
    void
    operator()(xmlpp::Node* pDimerizationGenNode) const throw(mzr::mzrXcpt)
    {
      // Cast the node pointer to element pointer; probably unnecessarily
      // dynamically.
      xmlpp::Element* pDimerizationGenElt
	= domUtils::mustBeElementPtr(pDimerizationGenNode);

      // Get the elements giving the binding partners.
      xmlpp::Node::NodeList molRefs
	= pDimerizationGenElt->get_children(plx::eltName::molRef);

      // Make sure there are two of them.
      if(2 != molRefs.size())
	throw domUtils::badChildCountXcpt(pDimerizationGenElt,
					  plx::eltName::molRef,
					  2,
					  (int) molRefs.size());

      // Parse the binding partners.
      std::vector<bindingPartner> bindingPartners(2);
      std::transform(molRefs.begin(),
		     molRefs.end(),
		     bindingPartners.begin(),
		     parseBindingPartner(rMolUnit));

      // Extract left binding partner's info.
      bnd::mol* pLeftMol = bindingPartners[0].first;
      int leftSiteNdx = bindingPartners[0].second;
      const bnd::bindingSite& rLeftSite = pLeftMol->getSiteByNdx(leftSiteNdx);
      const std::map<std::string, bnd::siteShape>& rLeftSiteShapeMap
	= rLeftSite.getSiteShapeMap();
      double leftMolWeight = pLeftMol->getDefaultParam()->getMolWeight();

      // Extract right binding partner's info.
      bnd::mol* pRightMol = bindingPartners[1].first;
      int rightSiteNdx = bindingPartners[1].second;
      const bnd::bindingSite& rRightSite = pRightMol->getSiteByNdx(rightSiteNdx);
      const std::map<std::string, bnd::siteShape>& rRightSiteShapeMap
	= rRightSite.getSiteShapeMap();
      double rightMolWeight = pRightMol->getDefaultParam()->getMolWeight();

      // Install the binding feature for the two binding sites.
      plx::bindingFeature* pBindingFeature
	= rPlexUnit.addBindingFeature(pLeftMol,
				       leftSiteNdx,
				       pRightMol,
				       rightSiteNdx);

      // Check if the user specified a reaction rate extrapolator.
      xmlpp::Attribute* pRateExtrapolatorAttr
	= pDimerizationGenElt->get_attribute
	(eltName::dimerizationGen_rateExtrapAttr);

      // Construct the dimerization reaction rate extrapolator according to
      // the option.  Reaction rate extrapolators are memory-managed by the
      // reaction generators they are attached to.
      dimerizeExtrapolator* pDimerizeExtrap = 0;
      if(0 != pRateExtrapolatorAttr
	 && (pRateExtrapolatorAttr->get_value() 
	     == eltName::dimerizationGen_rateExtrap_none))
	{
	  pDimerizeExtrap = new dimerizeNoExtrap();
	}
      else
	{
	  // This is for the "mass" case, for now the only alternative, as
	  // well as the default.
	  pDimerizeExtrap = new dimerizeMassExtrap(leftMolWeight,
						   rightMolWeight);
	}

      // There is only one option for decomposition rate extrapolation.
      decomposeExtrapolator* pDecompExtrap
	= new decomposeNoExtrap();

      // Parse the default rates.
      xmlpp::Element* pOnRateElt
	= domUtils::mustGetUniqueChild(pDimerizationGenElt,
				       eltName::defaultOnRate);
      double onRate
	= domUtils::mustGetAttrDouble(pOnRateElt,
				      eltName::defaultOnRate_valueAttr);
      xmlpp::Element* pOffRateElt
	= domUtils::mustGetUniqueChild(pDimerizationGenElt,
				       eltName::defaultOffRate);
      double offRate
	= domUtils::mustGetAttrDouble(pOffRateElt,
				      eltName::defaultOffRate_valueAttr);

      // Install the default rates as kinetics for all possible
      // combinations of site shapes for the two binding partners.
      //
      // Passing along the mols and binding sites for debugging.
      // NOTE: MOST PARAMETERS HERE ARE FOR DEBUGGING.
      std::for_each(rLeftSiteShapeMap.begin(),
		    rLeftSiteShapeMap.end(),
		    installDefaultKineticsToLeft(rRightSiteShapeMap,
						 onRate,
						 offRate,
						 pLeftMol,
						 rLeftSite,
						 pRightMol,
						 rRightSite,
						 pDimerizeExtrap,
						 pDecompExtrap));

      // Parse the allosteric shape pairs; i.e. the pairs of binding
      // site shapes whose kinetic differ from the default.
      //
      // NOTE: Now I shouldn't need the weights?
      xmlpp::Node::NodeList alloRatesNodes
	= pDimerizationGenElt->get_children(eltName::alloRates);
      std::for_each(alloRatesNodes.begin(),
		    alloRatesNodes.end(),
		    insertAlloRates(rLeftSite,
				    leftMolWeight,
				    rRightSite,
				    rightMolWeight,
				    pDimerizeExtrap,
				    pDecompExtrap));

      // Construct the decomposition reaction family, and add it to
      // the BAD BAD BAD static vector for memory management.
      decompFam* pDecompFam = new decompFam(rMzrUnit,
					    rPlexUnit,
					    pDecompExtrap);
      rMzrUnit.addReactionFamily(pDecompFam);

      // Attach the decomposition reaction generator to the binding feature.
      pBindingFeature->addRxnGen(pDecompFam->getRxnGen());

      // Construct the dimerizaton reaction family, and add it to
      // the BAD BAD BAD static vector for memory management.
      plx::siteFeature& rLeftSiteFeature
	= pLeftMol->getSiteFeature(leftSiteNdx);
      plx::siteFeature& rRightSiteFeature
	= pRightMol->getSiteFeature(rightSiteNdx);
      dimerizeFam* pDimerizeFam = new dimerizeFam(rLeftSiteFeature,
						  rRightSiteFeature,
						  rMzrUnit,
						  rPlexUnit,
						  pDimerizeExtrap);
      rMzrUnit.addReactionFamily(pDimerizeFam);

      // Attach the dimerization reaction generator the site features.
      rLeftSiteFeature.addRxnGen(pDimerizeFam->getLeftRxnGen());
      rRightSiteFeature.addRxnGen(pDimerizeFam->getRightRxnGen());
    }
  };

  void
  dimerUnit::parseDomInput(xmlpp::Element* pRootElement,
			   xmlpp::Element* pModelElement,
			   xmlpp::Element* pStreamsElement,
			   xmlpp::Element* pEventsElement) throw(std::exception)
  {
    // Get the header node for all reaction generators.
    xmlpp::Element* pReactionGensElt
      = domUtils::mustGetUniqueChild(pModelElement,
				     mzr::eltName::reactionGens);

    // Get all the dimerization-gen elements.  Eventually, I will have to
    // deal with the many reaction generators for the scaffold kinase reactions
    // and other special reactions.
    xmlpp::Node::NodeList dimerizationGens
      = pReactionGensElt->get_children(eltName::dimerizationGen);

    // Make dimerization and decomposition generators.
    std::for_each(dimerizationGens.begin(),
		  dimerizationGens.end(),
		  addDimerizeDecompose(rMzrUnit,
				       rMolUnit,
				       *this,
				       rPlexUnit));
  }
}
