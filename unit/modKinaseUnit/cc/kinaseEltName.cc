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

#include "modKinase/kinaseEltName.hh"

namespace kinase
{
  namespace eltName
  {
    const std::string kinaseGen("kinase-gen");
    const std::string modMolRef("mod-mol-ref");
    const std::string modMolRef_nameAttr("name");
    const std::string kinaseGen_rateExtrapAttr("rate-extrapolator");
    const std::string kinaseGen_rateExtrap_none("none");
    const std::string kinaseGen_rateExtrap_mass("mass");
    const std::string kinaseModMolRef("kinase-mod-mol-ref");
    const std::string kinaseModMolRef_nameAttr("name");
    const std::string atpModSiteRef("atp-mod-site-ref");
    const std::string atpModSiteRef_nameAttr("name");
    const std::string substrateModMolRef("substrate-mod-mol-ref");
    const std::string substrateModMolRef_nameAttr("name");
    const std::string phosModSiteRef("phos-mod-site-ref");
    const std::string phosModSiteRef_nameAttr("name");
    const std::string phosphorylatedModRef("phosphorylated-mod-ref");
    const std::string phosphorylatedModRef_nameAttr("name");
    const std::string atpBoundModRef("atp-bound-mod-ref");
    const std::string atpBoundModRef_nameAttr("name");
    const std::string adpBoundModRef("adp-bound-mod-ref");
    const std::string adpBoundModRef_nameAttr("name");

    const std::string ptaseGen("ptase-gen");
    const std::string ptaseGen_rateExtrapAttr("rate-extrapolator");
    const std::string ptaseGen_rateExtrap_none("none");
    const std::string ptaseGen_rateExtrap_mass("mass");
    const std::string ptaseStochSpeciesRef("ptase-stoch-species-ref");
    const std::string ptaseStochSpeciesRef_nameAttr("name");
    const std::string phosphateSpeciesRef("phosphate-species-ref");
    const std::string phosphateSpeciesRef_nameAttr("name");

    const std::string nucleotideBindGen("nucleotide-bind-gen");
    const std::string nucleotideBoundModRef("nucleotide-bound-mod-ref");
    const std::string nucleotideBoundModRef_nameAttr("name");
    const std::string noneModRef("none-mod-ref");
    const std::string noneModRef_nameAttr("name");
    const std::string nucleotideSpeciesRef("nucleotide-species-ref");
    const std::string nucleotideSpeciesRef_nameAttr("name");
    const std::string nucleotideBindGen_rateExtrapAttr("rate-extrapolator");
    const std::string nucleotideBindGen_rateExtrap_none("none");
    const std::string nucleotideBindGen_rateExtrap_mass("mass");
  }
}
