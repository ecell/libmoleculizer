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

#ifndef KINASEELTNAME_H
#define KINASEELTNAME_H

#include <string>

namespace kinase
{
  namespace eltName
  {
    extern const std::string kinaseGen;
    extern const std::string modMolRef;
    extern const std::string modMolRef_nameAttr;
    extern const std::string kinaseGen_rateExtrapAttr;
    extern const std::string kinaseGen_rateExtrap_none;
    extern const std::string kinaseGen_rateExtrap_mass;
    extern const std::string kinaseModMolRef;
    extern const std::string kinaseModMolRef_nameAttr;
    extern const std::string atpModSiteRef;
    extern const std::string atpModSiteRef_nameAttr;
    extern const std::string substrateModMolRef;
    extern const std::string substrateModMolRef_nameAttr;
    extern const std::string phosModSiteRef;
    extern const std::string phosModSiteRef_nameAttr;
    extern const std::string phosphorylatedModRef;
    extern const std::string phosphorylatedModRef_nameAttr;
    extern const std::string atpBoundModRef;
    extern const std::string atpBoundModRef_nameAttr;
    extern const std::string adpBoundModRef;
    extern const std::string adpBoundModRef_nameAttr;

    extern const std::string ptaseGen;
    extern const std::string ptaseGen_rateExtrapAttr;
    extern const std::string ptaseGen_rateExtrap_none;
    extern const std::string ptaseGen_rateExtrap_mass;
    extern const std::string ptaseStochSpeciesRef;
    extern const std::string ptaseStochSpeciesRef_nameAttr;
    extern const std::string phosphateSpeciesRef;
    extern const std::string phosphateSpeciesRef_nameAttr;

    extern const std::string nucleotideBindGen;
    extern const std::string nucleotideBoundModRef;
    extern const std::string nucleotideBoundModRef_nameAttr;
    extern const std::string noneModRef;
    extern const std::string noneModRef_nameAttr;
    extern const std::string nucleotideSpeciesRef;
    extern const std::string nucleotideSpeciesRef_nameAttr;
    extern const std::string nucleotideBindGen_rateExtrapAttr;
    extern const std::string nucleotideBindGen_rateExtrap_none;
    extern const std::string nucleotideBindGen_rateExtrap_mass;
  }
}

#endif
