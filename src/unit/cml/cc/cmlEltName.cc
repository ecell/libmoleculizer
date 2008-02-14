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

namespace cml
{
  namespace eltName
  {
    const std::string modifications("modifications");
    const std::string modification("modification");
    const std::string modification_nameAttr("name");
    const std::string weightDelta("weight-delta");
    const std::string weightDelta_daltonsAttr("daltons");

    const std::string mols("mols");
    const std::string modMol("mod-mol");
    const std::string modMol_nameAttr("name");
    const std::string weight("weight");
    const std::string weight_daltonsAttr("daltons");
    const std::string bindingSite("binding-site");
    const std::string bindingSite_nameAttr("name");
    const std::string defaultShapeRef("default-shape-ref");
    const std::string defaultShapeRef_nameAttr("name");
    const std::string siteShape("site-shape");
    const std::string siteShape_nameAttr("name");
    const std::string modSite("mod-site");
    const std::string modSite_nameAttr("name");
    const std::string defaultModRef("default-mod-ref");
    const std::string defaultModRef_nameAttr("name");

    const std::string smallMol("small-mol");
    const std::string smallMol_nameAttr("name");

    const std::string allostericState("allosteric-state");
    const std::string modMap("mod-map");
    const std::string modSiteRef("mod-site-ref");
    const std::string modSiteRef_nameAttr("name");
    const std::string modRef("mod-ref");
    const std::string modRef_nameAttr("name");
    const std::string siteShapeMap("site-shape-map");
    const std::string bindingSiteRef("binding-site-ref");
    const std::string bindingSiteRef_nameAttr("name");
    const std::string siteShapeRef("site-shape-ref");
    const std::string siteShapeRef_nameAttr("name");
  }
}
