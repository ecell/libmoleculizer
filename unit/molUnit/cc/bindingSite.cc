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

#include "mol/molEltName.hh"
#include "mol/bindingSite.hh"

namespace bnd
{
  class parseShapeForBindingSite :
    public std::unary_function<xmlpp::Node*,
    std::pair<std::string, siteShape> >
  {
  public:
    std::pair<std::string, siteShape>
    operator()(xmlpp::Node* pSiteShapeNode) throw(std::exception)
    {
      siteShape theShape(pSiteShapeNode);

      return std::make_pair(theShape.name,
			    theShape);
    }
  };
  
  bindingSite::bindingSite(xmlpp::Node* pBindingSiteNode)
    throw(std::exception)
  {
    xmlpp::Element* pBindingSiteElt
      = domUtils::mustBeElementPtr(pBindingSiteNode);

    name
      = domUtils::mustGetAttrString(pBindingSiteElt,
				    eltName::bindingSite_nameAttr);

    xmlpp::Node::NodeList siteShapeNodes
      = pBindingSiteElt->get_children(eltName::siteShape);
    std::transform(siteShapeNodes.begin(),
		   siteShapeNodes.end(),
		   std::inserter(shapesByName,
				 shapesByName.begin()),
		   parseShapeForBindingSite());

    xmlpp::Element* pDefaultShapeRefElt
      = domUtils::mustGetUniqueChild(pBindingSiteElt,
				     eltName::defaultShapeRef);

    std::string defaultShapeName
      = domUtils::mustGetAttrString
      (pDefaultShapeRefElt,
       eltName::defaultShapeRef_nameAttr);

    setDefaultParamByName(defaultShapeName);
  }

  siteParam
  bindingSite::
  mustGetParam(xmlpp::Node* pRequestingNode,
	       const mol* pMol,
	       const std::string& rShapeName) const
    throw(unknownSiteShapeXcpt)
  {
    siteParam theParam = getParam(rShapeName);

    if(! theParam)
      throw unknownSiteShapeXcpt(pRequestingNode,
				 *this,
				 pMol,
				 rShapeName);
    return theParam;
  }

  class insertShapeElt :
    public std::unary_function<std::pair<const std::string, siteShape>, void>
  {
    xmlpp::Element* pBindingSiteElt;
  public:
    insertShapeElt(xmlpp::Element* pBindingSiteElement) :
      pBindingSiteElt(pBindingSiteElement)
    {}

    void
    operator()(const std::pair<const std::string, siteShape>& rEntry) const
      throw(std::exception)
    {
      xmlpp::Element* pSiteShapeElt
	= pBindingSiteElt->add_child(eltName::siteShape);

      pSiteShapeElt->set_attribute(eltName::siteShape_nameAttr,
				   rEntry.second.name);
    }
  };

  xmlpp::Element*
  bindingSite::insertElt(xmlpp::Element* pMolElt) const throw(std::exception)
  {
    xmlpp::Element* pBindingSiteElt
      = pMolElt->add_child(eltName::bindingSite);

    pBindingSiteElt->set_attribute(eltName::bindingSite_nameAttr,
				   getName());

    // Put in the name of this binding site's default shape.
    xmlpp::Element* pDefaultShapeRefElt
      = pBindingSiteElt->add_child(eltName::defaultShapeRef);

    pDefaultShapeRefElt->set_attribute(eltName::defaultShapeRef_nameAttr,
				       getDefaultParam()->name);

    // Put in all the binding sites shapes.
    std::for_each(shapesByName.begin(),
		  shapesByName.end(),
		  insertShapeElt(pBindingSiteElt));

    return pBindingSiteElt;
  }
}
