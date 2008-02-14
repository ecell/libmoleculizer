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

#include "cml/cmlEltName.hh"
#include "cml/parseModSite.hh"
#include "utl/badElementCastXcpt.hh"
#include "cml/unkModXcpt.hh"

namespace cml
{
  std::pair<std::string, const cpx::modification*>
  parseModSite::
  operator()(xmlpp::Node* pModSiteNode) const 
    throw(utl::xcpt)
  {
    // Make sure the node is an element, possibly unnecessarily.
    xmlpp::Element* pModSiteElt
      = dynamic_cast<xmlpp::Element*>(pModSiteNode);
    if(0 == pModSiteElt) throw utl::dom::badElementCastXcpt(pModSiteNode);

    // Get the mod site name.
    std::string modSiteName
      = utl::dom::mustGetAttrString(pModSiteElt,
				    eltName::modSite_nameAttr);

    // Get the default modification element.
    xmlpp::Element* pDefaultModRefElt
      = utl::dom::mustGetUniqueChild(pModSiteElt,
				     eltName::defaultModRef);

    // Get the default modification name.
    std::string defaultModName
      = utl::dom::mustGetAttrString(pDefaultModRefElt,
				    eltName::defaultModRef_nameAttr);

    // Look up the default modification.
    // 
    // Note that this is yet another BAD BAD BAD static catalog use.
    // 
    // Testing that the lookup succeeded indicates lack of trust in
    // validation.
    const cpx::modification* pDefaultMod
      = rCmlUnit.getMod(defaultModName);
    if(0 == pDefaultMod) throw unkModXcpt(defaultModName,
					  pDefaultModRefElt);

    return std::make_pair(modSiteName,
			  pDefaultMod);
  }
}
