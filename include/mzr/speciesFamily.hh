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

#ifndef SPECIESFAMILY_H
#define SPECIESFAMILY_H

#include <map>
#include "mzr/notifyingSpecies.hh"

namespace mzr
{
  // This differs from the old speciesFamily mainly in not having
  // dumpables built in.

  template<class speciesKey, class memberSpecies>
  class speciesFamily :
    public speciesNotificationTarget<memberSpecies>,
    public std::map<speciesKey, memberSpecies*>
  {
  public:
    typedef
    typename speciesFamily<speciesKey, memberSpecies>::iterator
    iterator;

    // Construct memberSpecies from key, but doesn't install
    // it.  This is where plexFamily's allostery function lives,
    // for example.
    virtual memberSpecies*
    makeMember(const speciesKey& rKey) = 0;

    // Retrieves existing memberSpecies for key, or constructs
    // and installs result of makeMemeber.  Reaction generators
    // use this to construct their product species, for example.
    memberSpecies*
    getMember(const speciesKey& rKey);

    // Note that a speciesFamily must implement notify(pMemberSpecies)
    // in its role as notification proxy for all of its member species,
    // to notify all the targets that are interested in traits common
    // to all the species in the family.
  };

  template<class speciesKey, class memberSpecies>
  memberSpecies*
  speciesFamily<speciesKey, memberSpecies>::getMember(const speciesKey& rKey)
  {
    // Attempt to insert a bogus pointer.
    std::pair<iterator, bool> insertResult
      = insert(std::make_pair(rKey, (memberSpecies*) 0));

    // If the insert succeeded, then there is no species for the
    // presented key, and we construct and install a new one.
    if(insertResult.second)
      {
	// Make the new species.
	memberSpecies* pNewSpecies = makeMember(rKey);

	// Install it in the map.
	insertResult.first->second = pNewSpecies;

	return pNewSpecies;
      }
    else
      {
	return insertResult.first->second;
      }
  }
}

#endif // SPECIESFAMILY_H
