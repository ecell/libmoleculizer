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

#include "nucEx/nucExEltName.hh"

namespace nucEx
{
  namespace eltName
  {
    const std::string nucleotideExchangeGen("nucleotide-exchange-gen");
    const std::string nucleotideExchangeGen_rateExtrapAttr("rate-extrapolator");
    const std::string nucleotideExchangeGen_rateExtrap_none("none");
    const std::string nucleotideExchangeGen_rateExtrap_mass("mass");
    const std::string nucleotideBoundModRef("nucleotide-bound-mod-ref");
    const std::string nucleotideBoundModRef_nameAttr("name");
    const std::string noneModRef("none-mod-ref");
    const std::string noneModRef_nameAttr("name");
    const std::string nucleotideSpeciesRef("nucleotide-species-ref");
    const std::string nucleotideSpeciesRef_nameAttr("name");
    const std::string enabledOnRate("enabled-on-rate");
    const std::string enabledOnRate_valueAttr("value");
    const std::string plainOnRate("plain-on-rate");
    const std::string plainOnRate_valueAttr("value");
    const std::string enabledOffRate("enabled-off-rate");
    const std::string enabledOffRate_valueAttr("value");
    const std::string plainOffRate("plain-off-rate");
    const std::string plainOffRate_valueAttr("value");

    const std::string autoHydrolysisGen("auto-hydrolysis-gen");
    const std::string autoHydrolysisGen_rateExtrapAttr("rate-extrapolator");
    const std::string autoHydrolysisGen_rateExtrap_none("none");
    const std::string autoHydrolysisGen_rateExtrap_mass("mass");
    const std::string heteroHydrolysisGen("hetero-hydrolysis-gen");
    const std::string heteroHydrolysisGen_rateExtrapAttr("rate-extrapolator");
    const std::string heteroHydrolysisGen_rateExtrap_none("none");
    const std::string heteroHydrolysisGen_rateExtrap_mass("mass");
    const std::string unhydrolysedBoundModRef("unhydrolysed-bound-mod-ref");
    const std::string unhydrolysedBoundModRef_nameAttr("name");
    const std::string hydrolysedBoundModRef("hydrolysed-bound-mod-ref");
    const std::string hydrolysedBoundModRef_nameAttr("name");
  }
}
