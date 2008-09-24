//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

#ifndef MOL_PARSEMODMAP_H
#define MOL_PARSEMODMAP_H

#include "utl/dom.hh"
#include "cpx/modification.hh"
#include "mol/molUnit.hh"
#include "mol/unkModXcpt.hh"

namespace bnd
{
// This class is almost exactly the same as parseModSite;
// pretty much only the element names change.
//
// This is for parseing an allosteric state of a mod-mol.
class parseModMap : public
std::unary_function<xmlpp::Node*,
std::pair<std::string, const cpx::modification*> >
{
molUnit& rMolUnit;

public:
parseModMap(molUnit& refMolUnit) :
rMolUnit(refMolUnit)
{}

std::pair<std::string, const cpx::modification*>
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
const cpx::modification* pMod = rMolUnit.getMod(modName);
if(0 == pMod) throw unkModXcpt(modName,
pModRefElt);

return std::make_pair(modSiteName,
pMod);
}
};
}

#endif // MOL_PARSEMODMAP_H
