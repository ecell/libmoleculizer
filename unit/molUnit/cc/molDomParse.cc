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
#include "mzr/moleculizer.hh"
#include "mzr/mzrUnit.hh"
#include "mol/molXcpt.hh"
#include "mol/molDomParse.hh"
#include "mol/molUnit.hh"
#include "mol/modMol.hh"
#include "mol/smallMol.hh"
#include "mol/modQuery.hh"
#include "mol/molEltName.hh"
#include "plex/plexUnit.hh"

namespace bnd
{
  modMol*
  mustBeModMolPtr(xmlpp::Node* pRequestingNode,
		  mol* pMol)
    throw(badModMolCastXcpt)
  {
    modMol* pModMol = dynamic_cast<modMol*>(pMol);
    if(0 == pModMol) throw badModMolCastXcpt(pRequestingNode,
					     pMol);
    return pModMol;
  }

  smallMol*
  mustBeSmallMolPtr(xmlpp::Node* pRequestingNode,
		    mol* pMol)
    throw(badSmallMolCastXcpt)
  {
    smallMol* pSmallMol = dynamic_cast<smallMol*>(pMol);
    if(0 == pSmallMol) throw badSmallMolCastXcpt(pRequestingNode,
						 pMol);
    return pSmallMol;
  }

  class installModification :
    public std::unary_function<xmlpp::Node*, void>
  {
    molUnit& rMolUnit;
    
  public:
    installModification(molUnit& refMolUnit) :
      rMolUnit(refMolUnit)
    {}
    
//     void
//     operator()(xmlpp::Node* pModNode) const throw(mzr::mzrXcpt)
//     {
//       // When could I make this a static cast?
//       xmlpp::Element* pModElt
// 	= dynamic_cast<xmlpp::Element*>(pModNode);
//       if(0 == pModElt) throw domUtils::badElementCastXcpt(pModNode);

//       std::string modificationName
// 	= domUtils::mustGetAttrString
// 	(pModElt,
// 	 eltName::modification_nameAttr);

//       xmlpp::Element* pWeightDeltaElt
// 	= domUtils::mustGetUniqueChild(pModElt,
// 				       eltName::weightDelta);

//       double weightDelta
// 	= domUtils::mustGetAttrDouble
// 	(pWeightDeltaElt,
// 	 eltName::weightDelta_daltonsAttr);

//       // Install the BAD BAD BAD static catalog of the
//       // modification class.  internMod is okay, but it should
//       // be a molUnit method and should intern into an
//       // ordinary member of molUnit.
//       //
//       // I must decide whether this code will trust in validation or not.  It
//       // certainly makes some stuff simpler if I do.
//       if(rMolUnit.internMod(modificationName,
// 			    weightDelta))
// 	throw duplicateModNameXcpt(pModElt,
// 				   modificationName);
//     }

    void
    operator()(xmlpp::Node* pModNode) const throw(std::exception)
    {
      rMolUnit.mustAddParsedMod(pModNode);
    }
  };

//   class parseSiteShape :
//     public std::unary_function<xmlpp::Node*,
//     std::pair<std::string, siteShape> >
//   {
//   public:
//     std::pair<std::string, siteShape>
//     operator()(xmlpp::Node* pSiteShapeNode) const throw(mzr::mzrXcpt)
//     {
//       // Make sure the node is an element (possibly unnecessarily).
//       xmlpp::Element* pSiteShapeElt
// 	= dynamic_cast<xmlpp::Element*>(pSiteShapeNode);
//       if(0 == pSiteShapeElt)
// 	throw domUtils::badElementCastXcpt(pSiteShapeNode);

//       // Get the site shape name.
//       std::string siteShapeName
// 	= domUtils::mustGetAttrString(pSiteShapeElt,
// 				      eltName::siteShape_nameAttr);

//       // Test site shape name for uniqueness?
//       // At some point, I will have to decide whether I trust validation
//       // or not.

//       // Construct a site shape catalog entry.
//       return make_pair(siteShapeName,
// 		       siteShape(siteShapeName));
//     }
//   };

//   class parseBindingSite :
//     public std::unary_function<xmlpp::Node*, bindingSite>
//   {
//   public:
//     bindingSite
//     operator()(xmlpp::Node* pBindingSiteNode) const throw(mzr::mzrXcpt)
//     {
//       // When could I make this a static cast?
//       xmlpp::Element* pBindingSiteElt
// 	= dynamic_cast<xmlpp::Element*>(pBindingSiteNode);
//       if(0 == pBindingSiteElt)
// 	throw domUtils::badElementCastXcpt(pBindingSiteNode);

//       // Get the binding site name.
//       std::string siteName
// 	= domUtils::mustGetAttrString(pBindingSiteElt,
// 				      eltName::bindingSite_nameAttr);

//       // Process the site shapes.
//       std::map<std::string, siteShape> nameToShape;
//       xmlpp::Node::NodeList siteShapeNodes
// 	= pBindingSiteElt->get_children(eltName::siteShape);
//       std::transform(siteShapeNodes.begin(),
// 		     siteShapeNodes.end(),
// 		     std::inserter(nameToShape,
// 				   nameToShape.begin()),
// 		     parseSiteShape());

//       // Get the default shape name (it looks like the binding site
//       // constructor does the lookup using the above map.)
//       xmlpp::Element* pDefaultShapeRefElt
// 	= domUtils::mustGetUniqueChild(pBindingSiteElt,
// 				       eltName::defaultShapeRef);

//       std::string defaultShapeName
// 	= domUtils::mustGetAttrString
// 	(pDefaultShapeRefElt,
// 	 eltName::defaultShapeRef_nameAttr);

//       // Construct binding site.
//       return bindingSite(siteName,
// 			 nameToShape,
// 			 defaultShapeName);
//     }
//   };

  class parseModSite : public
  std::unary_function<xmlpp::Node*,
		       std::pair<std::string, const modification*> >
  {
    molUnit& rMolUnit;
    
  public:
    parseModSite(molUnit& refMolUnit) :
      rMolUnit(refMolUnit)
    {}
    
    std::pair<std::string, const modification*>
    operator()(xmlpp::Node* pModSiteNode) const throw(mzr::mzrXcpt)
    {
      // Make sure the node is an element, possibly unnecessarily.
      xmlpp::Element* pModSiteElt
	= dynamic_cast<xmlpp::Element*>(pModSiteNode);
      if(0 == pModSiteElt) throw domUtils::badElementCastXcpt(pModSiteNode);

      // Get the mod site name.
      std::string modSiteName
	= domUtils::mustGetAttrString(pModSiteElt,
				      eltName::modSite_nameAttr);

      // Get the default modification element.
      xmlpp::Element* pDefaultModRefElt
	= domUtils::mustGetUniqueChild(pModSiteElt,
				       eltName::defaultModRef);

      // Get the default modification name.
      std::string defaultModName
	= domUtils::mustGetAttrString(pDefaultModRefElt,
				      eltName::defaultModRef_nameAttr);

      // Look up the default modification.
      // 
      // Note that this is yet another BAD BAD BAD static catalog use.
      // 
      // Testing that the lookup succeeded indicates lack of trust in
      // validation.
      const modification* pDefaultMod
	= rMolUnit.getMod(defaultModName);
      if(0 == pDefaultMod) throw unknownModXcpt(pDefaultModRefElt,
						defaultModName);

      return std::make_pair(modSiteName,
			    pDefaultMod);
    }
  };

  // This class is almost exactly the same as parseModSite above;
  // pretty much only the element names change.
  class parseModMap : public
  std::unary_function<xmlpp::Node*,
		      std::pair<std::string, const modification*> >
  {
    molUnit& rMolUnit;
    
  public:
    parseModMap(molUnit& refMolUnit) :
      rMolUnit(refMolUnit)
    {}
    
    std::pair<std::string, const modification*>
    operator()(xmlpp::Node* pModSiteRefNode) const throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pModSiteRefElt
	= domUtils::mustBeElementPtr(pModSiteRefNode);

      // Get the mod site name.
      std::string modSiteName
	= domUtils::mustGetAttrString(pModSiteRefElt,
				      eltName::modSiteRef_nameAttr);

      // Get the mod-ref element that tells what modification is at this
      // modification site.
      xmlpp::Element* pModRefElt
	= domUtils::mustGetUniqueChild(pModSiteRefElt,
				       eltName::modRef);

      // Get the modification name.
      std::string modName
	= domUtils::mustGetAttrString(pModRefElt,
				      eltName::modRef_nameAttr);

      // Look up the modification in the BAD BAD BAD static catalog.
      //
      // Once again, testing for success here indicates lack of trust
      // in validation.
      const modification* pMod = rMolUnit.getMod(modName);
      if(0 == pMod) throw unknownModXcpt(pModRefElt,
					 modName);

      return std::make_pair(modSiteName,
			    pMod);
    }
  };

  class processAllostericSite :
    public std::unary_function<xmlpp::Node*, void>
  {
    molUnit& rMolUnit;
    std::vector<siteParam>& rParams;
    modMol* pMol;
  public:
    processAllostericSite(molUnit& refMolUnit,
			  modMol* pModMol,
			  std::vector<siteParam>& rSiteParams) :
      rMolUnit(refMolUnit),
      rParams(rSiteParams),
      pMol(pModMol)
    {}

    void
    operator()(xmlpp::Node* pBindingSiteRefNode) const throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pBindingSiteRefElt
	= domUtils::mustBeElementPtr(pBindingSiteRefNode);

      // Get the name of the binding site.
      std::string bindingSiteName
	= domUtils::mustGetAttrString
	(pBindingSiteRefElt,
	 eltName::bindingSiteRef_nameAttr);

      // Convert the binding site's name into the binding site's index.
      //
      // We need both the index and the binding site itself.
      int bindingSiteNdx;
      if(! pMol->findSite(bindingSiteName,
			  bindingSiteNdx))
	throw unknownSiteXcpt(pBindingSiteRefElt,
			      bindingSiteName);

      // Get the siteShapeRef element, which gives the allosteric shape
      // of the binding site.
      xmlpp::Element* pSiteShapeRefElt
	= domUtils::mustGetUniqueChild(pBindingSiteRefElt,
				       eltName::siteShapeRef);

      // Get the name of the site shape.
      std::string siteShapeName
	= domUtils::mustGetAttrString
	(pSiteShapeRefElt,
	 eltName::siteShapeRef_nameAttr);

      // Get pointer to the site shape.
      // 
      // (We need the binding site itself in order to get the site
      // shape, given its name.)
      const bindingSite& rBindingSite
	= pMol->getSiteByNdx(bindingSiteNdx);
      siteParam pSiteShape
	= rBindingSite.getParam(siteShapeName);
      if(! pSiteShape)
	throw unknownSiteShapeXcpt(pBindingSiteRefNode,
				   rBindingSite,
				   pMol,
				   siteShapeName);

      // Replace the default site shape of the binding site
      // with the specified site shape in the vector of site shapes.
      rParams[bindingSiteNdx] = pSiteShape;
    }
  };

  class installModMolAlloState :
    public std::unary_function<xmlpp::Node*, void>
  {
    molUnit& rMolUnit;
    modMol* pMol;
  public:
    installModMolAlloState(molUnit& refMolUnit,
			   modMol* pModMol) :
      rMolUnit(refMolUnit),
      pMol(pModMol)
    {}

    void
    operator()(xmlpp::Node* pAlloStateNode) const throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pAlloStateElt
	= domUtils::mustBeElementPtr(pAlloStateNode);

      // Get the map from modificaton sites to modifications that
      // describes this allosteric state.  This operation is very
      // similar to parsing the description of the modification
      // sites and their default modifications.
      //
      // It looks like the mod-map element could be eliminated.
      xmlpp::Element* pModMapElt
	= domUtils::mustGetUniqueChild(pAlloStateElt,
				       eltName::modMap);
      xmlpp::Node::NodeList modSiteRefs
	= pModMapElt->get_children(eltName::modSiteRef);
      std::map<std::string, const modification*> modMap;
      transform(modSiteRefs.begin(),
		modSiteRefs.end(),
		std::inserter(modMap,
			      modMap.begin()),
		parseModMap(rMolUnit));

      // Intern the modification map to get a state of the modMol.
      const modMolState* pState
	= pMol->internModMap(modMap);

      // Get the site shape map, which gives the shapes of the
      // mol's binding sites when in the above modification state.
      //
      // It looks like the site-shape-map element could be eliminated.
      xmlpp::Element* pSiteShapeMapElt
	= domUtils::mustGetUniqueChild(pAlloStateElt,
				       eltName::siteShapeMap);
      xmlpp::Node::NodeList bindingSiteRefs
	= pSiteShapeMapElt->get_children(eltName::bindingSiteRef);

      // Get the vector of site shapes currently associated with this
      // mol state.  (In all liklihood, this vector is actually generated
      // and installed in the mol's "allostery database" by this call to
      // alloMol::allostery.
      std::vector<siteParam>& rSiteShapes
	= pMol->allostery(pState);

      // Process the bindingSiteRefs, making modifications to the
      // defaultShapes, producing the allosteric vector of shapes.
      // Since this vector of siteParams is "hot," doing the replacements
      // through this 
      std::for_each(bindingSiteRefs.begin(),
		    bindingSiteRefs.end(),
		    processAllostericSite(rMolUnit,
					  pMol,
					  rSiteShapes));
    }
  };
    
  class installModMol :
    public std::unary_function<xmlpp::Node*, void>
  {

    molUnit& rMolUnit;
    
  public:
    installModMol(molUnit& refMolUnit) :
      rMolUnit(refMolUnit)
    {}
    
    void
    operator()(xmlpp::Node* pModMolNode) const throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pModMolElt
	= domUtils::mustBeElementPtr(pModMolNode);

      std::string modMolName
	= domUtils::mustGetAttrString(pModMolElt,
				      eltName::modMol_nameAttr);

      // Get the weight element.
      xmlpp::Element* pWeightElt
	= domUtils::mustGetUniqueChild(pModMolElt,
				       eltName::weight);

      // Get the weight.
      double weight
	= domUtils::mustGetAttrPosDouble(pWeightElt,
					 eltName::weight_daltonsAttr);

      // Process the binding sites.
      xmlpp::Node::NodeList bindingSiteNodes
	= pModMolElt->get_children(eltName::bindingSite);
      std::vector<bindingSite> bindingSites;
      std::transform(bindingSiteNodes.begin(),
		     bindingSiteNodes.end(),
		     std::back_inserter(bindingSites),
		     bindingSite::constructor());

      // Process the modification sites.
      xmlpp::Node::NodeList modSiteNodes
	= pModMolElt->get_children(eltName::modSite);
      std::map<std::string, const modification*> modMap;
      std::transform(modSiteNodes.begin(),
		     modSiteNodes.end(),
		     std::inserter(modMap,
				   modMap.begin()),
		     parseModSite(rMolUnit));

      // Construct the mol.  This has to be done before processing allosteric
      // states, because the mol's "allostery database" has to exist.
      modMol* pModMol
	= new modMol(modMolName,
		     bindingSites,
		     weight,
		     modMap);

      // Install the mol in plexUnit's catalog.  This would be okay if
      // plexUnit were allocated per-process, but as it is, this is
      // also BAD BAD BAD.
      rMolUnit.addMol(modMolName,
		      pModMol);

      // Process the allosteric states.
      xmlpp::Node::NodeList alloStateNodes
	= pModMolElt->get_children(eltName::allostericState);
      std::for_each(alloStateNodes.begin(),
		    alloStateNodes.end(),
		    installModMolAlloState(rMolUnit,
					   pModMol));
    }
  };

  class installSmallMol :
    public std::unary_function<xmlpp::Node*, void>
  {
    molUnit& rMolUnit;
  public:
    installSmallMol(molUnit& refMolUnit) :
      rMolUnit(refMolUnit)
    {
    }

    void
    operator()(xmlpp::Node* pSmallMolNode) const throw(mzr::mzrXcpt)
    {
      xmlpp::Element* pSmallMolElt
	= domUtils::mustBeElementPtr(pSmallMolNode);

      std::string smallMolName
	= domUtils::mustGetAttrString(pSmallMolElt,
				      eltName::smallMol_nameAttr);

      // Get the weight element.
      xmlpp::Element* pWeightElt
	= domUtils::mustGetUniqueChild(pSmallMolElt,
				       eltName::weight);

      // Get the molecular weight.
      double weight
	= domUtils::mustGetAttrPosDouble(pWeightElt,
					 eltName::weight_daltonsAttr);

      // Construct the smallMol.
      smallMol* pSmallMol
	= new smallMol(smallMolName,
		       weight);

      // Install the new smallMol in the catalog.
      rMolUnit.addMol(smallMolName,
		      pSmallMol);
    }
  };

  void
  molUnit::parseDomInput(xmlpp::Element* pRootElement,
			 xmlpp::Element* pModelElement,
			 xmlpp::Element* pStreamsElement,
			 xmlpp::Element* pEventsElement) throw(std::exception)
  {
    // Get the modifications element.
    xmlpp::Element* pModifications
      = domUtils::mustGetUniqueChild(pModelElement,
				     eltName::modifications);

    // Did some tests on the modification nodes to answer the following
    // question: do attributes show up in get_children results?
    //
    // This seems possible, since a list of arbitrary Node*'s is returned.  But
    // only elements result from Node::add_child, so libxml++ may regard
    // Elements as the only legitimate "children."
    //
    // So far, it looks like the only reason for not returning a list of
    // elements is that they don't want to have a new type for it.  NodeList is
    // also used to return results from xpath queries, and therefore must be
    // able to contain Node::Attribute*'s.

    // Install the modifications.  For now, these are going into a static
    // catalog in the modMixin class.  (Note that "installModification" does not
    // require "pTheMolculizer" as an argument, a very very bad sign.)  This
    // catalog should be an ordinary member of the mol unit, and units should be
    // allocated by new DURING MODULE LOADING, rather than in .so
    // initialization, so that different moleculizer processes will get
    // different module objects.

    xmlpp::Node::NodeList mods
      = pModifications->get_children(eltName::modification);
    std::for_each(mods.begin(),
		  mods.end(),
		  installModification(*this));

    // Get the mols element.
    xmlpp::Element* pMols
      = domUtils::mustGetUniqueChild(pModelElement,
				     eltName::mols);

    // Install the small-mols, which represent small molecules that
    // participate in binding (e.g. ATP) and therefore need the binding
    // machinery that mols have.
    xmlpp::Node::NodeList smallMols
      = pMols->get_children(eltName::smallMol);
    std::for_each(smallMols.begin(),
		  smallMols.end(),
		  installSmallMol(*this));

    // Install the mod-mols.  Again note that, due to the evil static
    // catalogs, etc. the argument list to "installModMol" is way too short.
    xmlpp::Node::NodeList modMols
      = pMols->get_children(eltName::modMol);
    std::for_each(modMols.begin(),
		  modMols.end(),
		  installModMol(*this));
  }
}
