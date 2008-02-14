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

#include "clx/clxUnit.hh"
#include "clx/dupNodeOmniXcpt.hh"
#include "clx/noOmniForNodeXcpt.hh"

namespace clx
{
  clxUnit::
  clxUnit(cpt::cptApp& rCptApp,
	  cpt::cptUnit& refCptUnit,
	  cml::cmlUnit& refCmlUnit) :
    cpt::unit("plex",
	      rCptApp),
    rCptUnit(refCptUnit),
    rCmlUnit(refCmlUnit),
    recognize(*this,
	      refCptUnit,
	      refCptUnit.getCompartmentGraph())
  {
    // Model elements for which plex unit is responsible.
    inputCap.addModelContentName(eltName::allostericPlexes);
    inputCap.addModelContentName(eltName::allostericOmnis);

    // Explicit species for which plex unit is responsible.
    inputCap.addExplicitSpeciesContentName(eltName::plexSpecies);

    // Species streams for which plex unit is responsible.
    inputCap.addSpeciesStreamsContentName
      (eltName::plexSpeciesStream);
    inputCap.addSpeciesStreamsContentName
      (eltName::omniSpeciesStream);
    
    // Not responsible for any reaction generators.
    // Not responsible for any species.

    // The clxUnit may as well register its own omniplex Xpaths.
    const std::string slash("/");

    std::ostringstream alloOmnisXpath;
    alloOmnisXpath << cpt::eltName::model
		   << slash
		   << eltName::allostericOmnis
		   << slash
		   << eltName::allostericOmni;
    addOmniXpath(alloOmnisXpath.str());

    std::ostringstream omniSpeciesStreamXpath;
    omniSpeciesStreamXpath << cpt::eltName::streams
			   << slash
			   << cpt::eltName::speciesStreams
			   << slash
			   << eltName::omniSpeciesStream;
    addOmniXpath(omniSpeciesStreamXpath.str());
  }

  // Installs a new binding feature for the given pair of structural
  // sites. Any plex that binds the two given structural bindings
  // together will need use the bindingFeature associated to the pair
  // of structural sites to inform the family of decompositions of
  // the structural binding.
  fnd::feature<cpx::cxBinding<cptPlexSpecies, cptPlexFamily> >*
  clxUnit::addBindingFeature(cml::cptMol* pLeftMol,
			     int leftMolSiteSpec,
			     cml::cptMol* pRightMol,
			     int rightMolSiteSpec)
  {
    cpx::structuralSite<cml::cptMol> leftSite(pLeftMol,
					      leftMolSiteSpec);

    cpx::structuralSite<cml::cptMol> rightSite(pRightMol,
					       rightMolSiteSpec);

    cpx::structuralBinding<cml::cptMol> sBinding(leftSite,
						 rightSite);

    fnd::feature<cpx::cxBinding<cptPlexSpecies, cptPlexFamily> > 
      bindingFtr;

    // I note (13Jul05) that I don't take exception to failure of insertion.
    cpx::knownBindings<cml::cptMol, fnd::feature<cpx::cxBinding<cptPlexSpecies, cptPlexFamily> > >::iterator
      iEntry = bindingFeatures.insert(std::make_pair(sBinding,
						     bindingFtr)).first;

    return &(iEntry->second);
  }

  // I expect this routine to be used in the constructor for plexFamily,
  // to construct the plexFamily's feature map for bindings.
  //
  // This routine looks for the feature "in both directions."
  fnd::feature<cpx::cxBinding<cptPlexSpecies, cptPlexFamily> >*
  clxUnit::findBindingFeature(cml::cptMol* pLeftMol,
			      int leftMolSiteSpec,
			      cml::cptMol* pRightMol,
			      int rightMolSiteSpec)
  {
    cpx::structuralSite<cml::cptMol> leftSite(pLeftMol,
					      leftMolSiteSpec);

    cpx::structuralSite<cml::cptMol> rightSite(pRightMol,
					       rightMolSiteSpec);

    cpx::structuralBinding<cml::cptMol> sBinding(leftSite,
						 rightSite);

    cpx::knownBindings<cml::cptMol, fnd::feature<cpx::cxBinding<cptPlexSpecies, cptPlexFamily> > >::iterator
      iEntry = bindingFeatures.find(sBinding);

    // May have to try the sites in the opposite order.
    if(iEntry == bindingFeatures.end())
      {
	cpx::structuralBinding<cml::cptMol> oppBinding(rightSite,
						       leftSite);
	iEntry = bindingFeatures.find(oppBinding);
      }

    return iEntry == bindingFeatures.end()
      ? 0
      : &(iEntry->second);
  }

  void
  clxUnit::
  addOmniPlex(cptOmniPlex* pOmniPlex,
	      xmlpp::Node* pParentNode)
    throw(utl::xcpt)
  {
    // Many omniplexes may have the same plexFamily, so
    // we should pay no attention if this insertion fails.
    omniPlexFamilies.insert(pOmniPlex->getFamily());

    // If there is already an entry for the node, then something
    // is internally exceptional.
    if(! (nodeToOmni.insert(std::make_pair(pParentNode,
					   pOmniPlex)).second))
      throw dupNodeOmniXcpt(pParentNode);
  }

  cptOmniPlex*
  clxUnit::
  getOmniForNode(xmlpp::Node* pParentNode) const
  {
    std::map<xmlpp::Node*, cptOmniPlex*>::const_iterator
      iEntry
      = nodeToOmni.find(pParentNode);

    return (nodeToOmni.end() == iEntry)
      ? 0
      : iEntry->second;
  }

  cptOmniPlex*
  clxUnit::
  mustGetOmniForNode(xmlpp::Node* pParentNode) const
    throw(utl::xcpt)
  {
    cptOmniPlex* pOmni
      = getOmniForNode(pParentNode);

    if(! pOmni) throw noOmniForNodeXcpt(pParentNode);

    return pOmni;
  }

  // For dumping state in XML.
  void
  clxUnit::
  insertStateElts(xmlpp::Element* pRootElt)
    throw(std::exception)
  {
    // Get the model element.
    xmlpp::Element* pModelElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     cpt::eltName::model);

    // Ensure that the tagged-species element is here.
    xmlpp::Element* pTaggedSpeciesElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     cpt::eltName::taggedSpecies);

    // Insert tagged-plex-species nodes.
    recognize.insertSpecies(pTaggedSpeciesElt);
  }
}
