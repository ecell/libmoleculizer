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

#include "bndKinase/bndKinaseEltName.hh"

namespace bndKinase
{
  namespace eltName
  {
    // Kinase reactions with nucleotide binding treated as ordinary binding.
    const std::string bndKinaseGen("bnd-kinase-gen");
    const std::string substrateModMolInstanceRef("substrate-mod-mol-instance-ref");
    const std::string substrateModMolInstanceRef_nameAttr("name");
    const std::string phosSiteRef("phos-site-ref");
    const std::string phosSiteRef_nameAttr("name");
    const std::string phosphorylatedModRef("phosphorylated-mod-ref");
    const std::string phosphorylatedModRef_nameAttr("name");
    const std::string atpSmallMolInstanceRef("atp-small-mol-instance-ref");
    const std::string atpSmallMolInstanceRef_nameAttr("name");
    const std::string adpSmallMolRef("adp-small-mol-ref");
    const std::string adpSmallMolRef_nameAttr("name");
    const std::string rate("rate");
    const std::string rate_valueAttr("value");

    // Generic modification reaction generator.
    const std::string modGen("mod-gen");
    const std::string modSiteRef("mod-site-ref");
    const std::string modSiteRef_nameAttr("name");
    const std::string installedModRef("installed-mod-ref");
    const std::string installedModRef_nameAttr("name");
  }
}
