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

#ifndef ODIEELTNAME_H
#define ODIEELTNAME_H

#include <string>

namespace odie
{
  namespace odieEltName
  {
    extern const std::string odie;
  
    extern const std::string species;
    extern const std::string speciesEntry;
    extern const std::string speciesEntry_nameAttr;
    extern const std::string speciesEntry_concAttr;

    extern const std::string reactions;
    extern const std::string reaction;
    extern const std::string reactionSubstrate;
    extern const std::string reactionSubstrate_nameAttr;
    extern const std::string reactionSubstrate_multAttr;
    extern const std::string reactionProduct;
    extern const std::string reactionProduct_nameAttr;
    extern const std::string reactionProduct_multAttr;
    extern const std::string rate;
    extern const std::string rate_valueAttr;

    extern const std::string dumpables;
    extern const std::string dumpable;
    extern const std::string dumpable_nameAttr;
    extern const std::string speciesRef;
    extern const std::string speciesRef_nameAttr;

    extern const std::string dumpStreams;
    extern const std::string dumpStream;
    extern const std::string dumpStream_nameAttr;
    extern const std::string dumpableRef;
    extern const std::string dumpableRef_nameAttr;

    extern const std::string stopTime;
    extern const std::string stopTime_secondsAttr;

    extern const std::string dumpInterval;
    extern const std::string dumpInterval_secondsAttr;

    extern const std::string epsilonAbs;
    extern const std::string epsilonAbs_valueAttr;
    extern const std::string epsilonRel;
    extern const std::string epsilonRel_valueAttr;
    extern const std::string stateScalar;
    extern const std::string stateScalar_valueAttr;
    extern const std::string derivScalar;
    extern const std::string derivScalar_valueAttr;
  };
}

#endif // ODIEELTNAME_H
