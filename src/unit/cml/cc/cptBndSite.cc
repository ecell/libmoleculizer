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

#include "cml/cptBndSite.hh"
#include "cml/cmlEltName.hh"
#include "cml/unkSiteShapeXcpt.hh"
#include "clx/cptPlexFamily.hh"

namespace cml
{
  cptBndSite::
  cptBndSite(const std::string& rName,
	     const std::set<std::string>& rShapeNames,
	     const std::string& rDefaultShapeName) :
    cpx::basicBndSite(rName,
		      rShapeNames,
		      rDefaultShapeName)
  {}

  const cpx::siteShape*
  cptBndSite::
  mustGetShape(const cptMol* pMol,
	       const std::string& rShapeName,
	       const xmlpp::Node* pRequestingNode) const
    throw(utl::xcpt)
  {
    const cpx::siteShape* pShape
      = getShape(rShapeName);

    if(! pShape)
      throw unkSiteShapeXcpt(pRequestingNode,
			     *this,
			     pMol,
			     rShapeName);

    return pShape;
  }

  class insertShapeElt :
    public std::unary_function<std::pair<const std::string, cpx::siteShape>, void>
  {
    xmlpp::Element* pBindingSiteElt;
  public:
    insertShapeElt(xmlpp::Element* pBindingSiteElement) :
      pBindingSiteElt(pBindingSiteElement)
    {}

    void
    operator()(const std::pair<const std::string, cpx::siteShape>& rEntry) const
      throw(std::exception)
    {
      xmlpp::Element* pSiteShapeElt
	= pBindingSiteElt->add_child(eltName::siteShape);

      pSiteShapeElt->set_attribute(eltName::siteShape_nameAttr,
				   rEntry.second.getName());
    }
  };

  xmlpp::Element*
  cptBndSite::
  insertElt(xmlpp::Element* pMolElt) const 
    throw(utl::xcpt)
  {
    xmlpp::Element* pBindingSiteElt
      = pMolElt->add_child(eltName::bindingSite);

    pBindingSiteElt->set_attribute(eltName::bindingSite_nameAttr,
				   getName());

    // Put in the name of this binding site's default shape.
    xmlpp::Element* pDefaultShapeRefElt
      = pBindingSiteElt->add_child(eltName::defaultShapeRef);

    pDefaultShapeRefElt->set_attribute(eltName::defaultShapeRef_nameAttr,
				       getDefaultShape()->getName());

    // Put in all the binding sites shapes.
    std::for_each(shapesByName.begin(),
		  shapesByName.end(),
		  insertShapeElt(pBindingSiteElt));

    return pBindingSiteElt;
  }
}
