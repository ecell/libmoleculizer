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

#ifndef BNDKINASEELTNAME_H
#define BNDKINASEELTNAME_H

#include <string>

namespace bndKinase
{
  namespace eltName
  {
    // Kinase reactions with nucleotide binding treated as ordinary binding.
    extern const std::string bndKinaseGen;
    extern const std::string substrateModMolInstanceRef;
    extern const std::string substrateModMolInstanceRef_nameAttr;
    extern const std::string phosSiteRef;
    extern const std::string phosSiteRef_nameAttr;
    extern const std::string phosphorylatedModRef;
    extern const std::string phosphorylatedModRef_nameAttr;
    extern const std::string atpSmallMolInstanceRef;
    extern const std::string atpSmallMolInstanceRef_nameAttr;
    extern const std::string adpSmallMolRef;
    extern const std::string adpSmallMolRef_nameAttr;
    extern const std::string rate;
    extern const std::string rate_valueAttr;

    // Generic modification reaction generator.
    extern const std::string modGen;
    extern const std::string modSiteRef;
    extern const std::string modSiteRef_nameAttr;
    extern const std::string installedModRef;
    extern const std::string installedModRef_nameAttr;
  }
}

#endif // BNDKINASEELTNAME_H
