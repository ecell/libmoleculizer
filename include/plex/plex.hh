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

#ifndef PLEX_H
#define PLEX_H

/*! \defgroup plexSpeciesGroup Complexes
  \ingroup plexGroup

  \brief Protein complexes and their constituents. */

#include <vector>
#include <map>
#include <string>
#include "mol/mol.hh"
#include "plex/plexSpec.hh"
#include "plex/plexXcpt.hh"

namespace plx
{
  class plexIsoPair;

  /*! \ingroup plexSpeciesGroup
    \brief A protein complex. */
  class plex
  {
    class insertBindingElt;
  
  public:
    std::vector<bnd::mol*> mols;		// essentially graph vertices.
    std::vector<plexBinding> bindings;	// essentially graph edges.

    // Used in finding the free sites on a plex.
    void
    makeSiteToBindings(std::map<plexSiteSpec, int>& rSiteToBindings) const;

    // Finds the free sites on a plex.
    void
    makeFreeSiteVector(std::vector<plexSiteSpec>& rFreeSiteVector) const;

    // Used in finding connected components of a plex.
    void
    makeMolToBindings(std::multimap<int, int>& rMolToBindings) const;

    // Used in "topological hashing" of plexes.
    void
    makeMolToMol(std::multimap<bnd::mol*, bnd::mol*>& rMolToMol) const;
  
    // Used in finding connected components of a plex.
    void
    pushConnectedBindings(int molNdx,
			  const std::multimap<int, int>& rMolToBindings,
			  plexIsoPair& rIso,
			  plex& component) const;
  
    // Used in finding connected components of a plex.
    int
    pushConnectedMol(int molNdx,
		     plexIsoPair& rIso,
		     plex& component) const;

    // Loads the connected component of molNdx into component, which should
    // be empty when the routine is called.  An partial isomoprhism from the
    // subcomplex of this plex onto the component is recorded in rIso. That
    // plexIsoPair should be large enough to hold the entire plex, in case
    // this plex is connected.  Create this isomorphism as is done in
    // makeConnectedComponent:
    // plexIsoPair iso(mols.size(), bindings.size());
    void
    makeTrackedComponent(int molNdx,
			 plex& component,
			 plexIsoPair& rIso) const;

    // This does the same as above, but creates and disposes of an appropriate
    // sized isomorphism.
    void
    makeConnectedComponent(int molNdx,
			   plex& component) const;

    int
    hashValue(void) const;

    // Used in the recognition cache.
    bool
    operator<(const plex& rRightPlex) const;

    // Plexes occur in various contexts.  For now, I intend to just
    // use stringified instance indices as instance names.
    xmlpp::Element*
    insertElt(xmlpp::Element* pParentElt) const
      throw(std::exception);
  };
}

#endif
