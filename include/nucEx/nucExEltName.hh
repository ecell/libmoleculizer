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

#ifndef NUCEXELTNAME_H
#define NUCEXELTNAME_H

#include <string>

namespace nucEx
{
  namespace eltName
  {
    extern const std::string nucleotideExchangeGen;
    extern const std::string nucleotideExchangeGen_rateExtrapAttr;
    extern const std::string nucleotideExchangeGen_rateExtrap_none;
    extern const std::string nucleotideExchangeGen_rateExtrap_mass;
    extern const std::string nucleotideBoundModRef;
    extern const std::string nucleotideBoundModRef_nameAttr;
    extern const std::string noneModRef;
    extern const std::string noneModRef_nameAttr;
    extern const std::string nucleotideSpeciesRef;
    extern const std::string nucleotideSpeciesRef_nameAttr;
    extern const std::string enabledOnRate;
    extern const std::string enabledOnRate_valueAttr;
    extern const std::string plainOnRate;
    extern const std::string plainOnRate_valueAttr;
    extern const std::string enabledOffRate;
    extern const std::string enabledOffRate_valueAttr;
    extern const std::string plainOffRate;
    extern const std::string plainOffRate_valueAttr;

    extern const std::string autoHydrolysisGen;
    extern const std::string autoHydrolysisGen_rateExtrapAttr;
    extern const std::string autoHydrolysisGen_rateExtrap_none;
    extern const std::string autoHydrolysisGen_rateExtrap_mass;
    extern const std::string heteroHydrolysisGen;
    extern const std::string heteroHydrolysisGen_rateExtrapAttr;
    extern const std::string heteroHydrolysisGen_rateExtrap_none;
    extern const std::string heteroHydrolysisGen_rateExtrap_mass;
    extern const std::string unhydrolysedBoundModRef;
    extern const std::string unhydrolysedBoundModRef_nameAttr;
    extern const std::string hydrolysedBoundModRef;
    extern const std::string hydrolysedBoundModRef_nameAttr;
  }
}

#endif
