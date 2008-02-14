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
#include "cml/parseModMap.hh"
#include "cml/unkModXcpt.hh"

namespace cml
{
  std::pair<std::string, const cpx::modification*>
  parseModMap::
  operator()(xmlpp::Node* pModSiteRefNode) const 
    throw(utl::xcpt)
  {
    xmlpp::Element* pModSiteRefElt
      = utl::dom::mustBeElementPtr(pModSiteRefNode);

    // Get the mod site name.
    std::string modSiteName
      = utl::dom::mustGetAttrString(pModSiteRefElt,
				    eltName::modSiteRef_nameAttr);

    // Get the mod-ref element that tells what modification is at this
    // modification site.
    xmlpp::Element* pModRefElt
      = utl::dom::mustGetUniqueChild(pModSiteRefElt,
				     eltName::modRef);

    // Get the modification name.
    std::string modName
      = utl::dom::mustGetAttrString(pModRefElt,
				    eltName::modRef_nameAttr);

    // Look up the modification in the BAD BAD BAD static catalog.
    //
    // Once again, testing for success here indicates lack of trust
    // in validation.
    const cpx::modification* pMod = rCmlUnit.getMod(modName);
    if(0 == pMod) throw unkModXcpt(modName,
				   pModRefElt);

    return std::make_pair(modSiteName,
			  pMod);
  }
}
