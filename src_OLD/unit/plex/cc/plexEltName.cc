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

#include "plex/plexEltName.hh"

namespace plx
{
  namespace eltName
  {
    const std::string allostericPlexes("allosteric-plexes");
    const std::string allostericPlex("allosteric-plex");
    const std::string allostericOmnis("allosteric-omnis");
    const std::string allostericOmni("allosteric-omni");
    const std::string plex("plex");
    const std::string molInstance("mol-instance");
    const std::string molInstance_nameAttr("name");
    const std::string molRef("mol-ref");
    const std::string molRef_nameAttr("name");
    const std::string binding("binding");
    const std::string molInstanceRef("mol-instance-ref");
    const std::string molInstanceRef_nameAttr("name");
    const std::string allostericSites("allosteric-sites");

    const std::string plexSpecies("plex-species");
    const std::string plexSpecies_nameAttr("name");

    const std::string omniSpeciesStream("omni-species-stream");
    const std::string omniSpeciesStream_nameAttr("name");
    const std::string plexSpeciesStream("plex-species-stream");
    const std::string plexSpeciesStream_nameAttr("name");

    const std::string instanceStates("instance-states");
    const std::string modMolInstanceRef("mod-mol-instance-ref");
    const std::string modMolInstanceRef_nameAttr("name");

    const std::string unboundSites("unbound-sites");
    const std::string instanceRef("instance-ref");
    const std::string instanceRef_nameAttr("name");
    const std::string siteRef("site-ref");
    const std::string siteRef_nameAttr("name");

    // For state dump.

    const std::string taggedPlexSpecies("tagged-plex-species");
    const std::string taggedPlexSpecies_tagAttr("tag");
    const std::string taggedPlexSpecies_nameAttr("name");
    const std::string population("population");
    const std::string population_countAttr("count");
    const std::string concentration("concentration");
    const std::string concentration_valueAttr("value");
    const std::string updated("updated");
  }
}
