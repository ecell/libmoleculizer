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

#include "rk4/tauEltName.hh"

namespace rk4tau
{
  namespace tauEltName
  {
    const std::string rk4tau("rk4tau");

    const std::string species("species");
    const std::string speciesEntry("species-entry");
    const std::string speciesEntry_nameAttr("name");
    const std::string speciesEntry_popAttr("pop");
  
    const std::string reactions("reactions");
    const std::string reaction("reaction");
    const std::string reactionSubstrate("reaction-substrate");
    const std::string reactionSubstrate_nameAttr("name");
    const std::string reactionSubstrate_multAttr("multiplicity");
    const std::string reactionProduct("reaction-product");
    const std::string reactionProduct_nameAttr("name");
    const std::string reactionProduct_multAttr("multiplicity");
    const std::string rate("rate");
    const std::string rate_valueAttr("value");

    const std::string dumpables("dumpables");
    const std::string dumpable("dumpable");
    const std::string dumpable_nameAttr("name");
    const std::string speciesRef("species-ref");
    const std::string speciesRef_nameAttr("name");

    const std::string dumpStreams("dump-streams");
    const std::string dumpStream("dump-stream");
    // Should this be "fileName" ?
    const std::string dumpStream_nameAttr("name");
    const std::string dumpableRef("dumpable-ref");
    const std::string dumpableRef_nameAttr("name");

    const std::string volume("volume");
    const std::string volume_litersAttr("liters");
    const std::string stopTime("stop-time");
    const std::string stopTime_secondsAttr("seconds");
    const std::string dumpInterval("dump-interval");
    const std::string dumpInterval_secondsAttr("seconds");
    const std::string epsilon("epsilon");
    const std::string epsilon_valueAttr("value");
    const std::string seed("seed");
    const std::string seed_valueAttr("value");
  }
}
