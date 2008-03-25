/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008  Walter Lawrence (Larry) Lok.
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

#ifndef MOL_PARSEBNDSITE_H
#define MOL_PARSEBNDSITE_H

#include "utl/dom.hh"
#include "mol/mzrBndSite.hh"
#include "mol/parseSiteShapeName.hh"
#include "mol/molEltName.hh"

namespace bnd
{
  class parseMzrBndSite :
    public std::unary_function<xmlpp::Node*, mzrBndSite>
  {
  public:
    mzrBndSite
    operator()(xmlpp::Node* pBindingSiteNode) const
      throw(utl::xcpt)
    {
      xmlpp::Element* pBindingSiteElt
	= utl::dom::mustBeElementPtr(pBindingSiteNode);

      std::string name
	= utl::dom::mustGetAttrString(pBindingSiteElt,
				      eltName::bindingSite_nameAttr);


      std::set<std::string> siteShapeNames;
      xmlpp::Node::NodeList siteShapeNodes
	= pBindingSiteElt->get_children(eltName::siteShape);
      std::transform(siteShapeNodes.begin(),
		     siteShapeNodes.end(),
		     std::inserter(siteShapeNames,
				   siteShapeNames.begin()),
		     parseSiteShapeName());

      xmlpp::Element* pDefaultShapeRefElt
	= utl::dom::mustGetUniqueChild(pBindingSiteElt,
				       eltName::defaultShapeRef);

      std::string defaultShapeName
	= utl::dom::mustGetAttrString(pDefaultShapeRefElt,
				      eltName::defaultShapeRef_nameAttr);

      return mzrBndSite(name,
			siteShapeNames,
			defaultShapeName);
    }
  };
}

#endif // MOL_PARSEBNDSITE_H
