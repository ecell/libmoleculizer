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

#ifndef RECOGNIZER_H
#define RECOGNIZER_H

#include <map>
#include <vector>
#include <fstream>
#include "plex/plexMap.hh"

namespace plx
{
  class plexFamily;
  class plexUnit;

  /*! \ingroup plexStructGroup
    \brief Class to maintain structural classification of complexes. */
  class recognizer
  {
    // In some cases, we have to know the isomorphism between
    // the recognized plex and the paradigm of its family.  This
    // is computed when the plex is first recognized, and saved
    // in the cache.
    class recognition
    {
    public:
      plexFamily* pPlexFamily;
      plexIsoPair iso;
    };
    std::map<plex, recognition> recognizedCache;

    plexUnit& rPlexUnit;

  public:

    recognizer(plexUnit& refPlexUnit) :
      rPlexUnit(refPlexUnit)
    {}

    // Publicized in order to traverse all the plexFamilies.
    //
    // In particular, for plexUnit::prepareToRun().
    std::multimap<int, plexFamily*> plexHasher;

    ~recognizer(void);

    int familyCount(void) const
    {
      return plexHasher.size();
    }

    // Finds the plexFamily of a plex, and gives the isomorphism of the
    // given plex with the plexFamily's paradigm.  This just runs the
    // bare constructor of the plexFamily, leaving undone the "phase
    // 2" part of plexFamily initialization: connection to features and
    // generation of the default parameter.
    //
    // I expect this to be used during setup, rather than at runtime.
    // Its purpose is to give finer control over the construction of
    // the plexFamily for a user-defined complex, as opposed to an
    // automatically generated one.  The user can make arbitrary
    // allosteric modifications, so we need (?) a way of unifying
    // families before constructing their default parameter/species.

    // This is a replacement function for the above.  It has an unpleasant
    // interface, but it thereby avoids some replication of code.
    bool
    unify(const plex& rPlex,
	  plexFamily*& rpFamily,
	  plexIsoPair* pIso = 0);

    // Recognizes an ordinary plex, and produces a fully initialized
    // plexFamily.  I expect this to be used at runtime.
    plexFamily*
    operator()(const plex& aPlex);

    // Recongizes an ordinary plex, and produces a fully initialized
    // plexFamily, as the above, but also reports the isomorphism from
    // the given plex to the plexFamily's paradigm.
    //
    // This is used in the decomposition reaction.
    plexFamily*
    operator()(const plex& aPlex,
	       plexIsoPair& rIso);

    // Output routine.
    void
    insertSpecies(xmlpp::Element* pExplicitSpeciesElt,
		  double molarFactor) const
      throw(std::exception);
  };
}

#endif
