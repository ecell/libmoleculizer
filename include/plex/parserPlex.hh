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

#ifndef PARSERPLEX_H
#define PARSERPLEX_H

#include "plex/plex.hh"

namespace plx
{
    /*! \ingroup plexSpeciesGroup
    \brief Complex embellished with instance names for the mols.

    This is for use in the parser: you need to be able to
    refer to the constituent mols in some easier way than, say,
    indices. */
  class parserPlex : public plex
  {
  public:
    // For incremental construction.
    parserPlex(void)
    {}

    // For constructing a parser plex (which has instance names) from
    // a plex and a naming of the instances.
    parserPlex(const std::map<std::string, int>& rNameToMolNdx,
	       const plex& rOriginal) :
      plex(rOriginal),
      nameToMolNdx(rNameToMolNdx)
    {}

    // Mapping of instance names to mol index.
    std::map<std::string, int> nameToMolNdx;

    // Appends a mol with the given instance name, returning
    // the new instance's index.
    int
    addMolByName(const std::string& rName,
		 bnd::mol* pMol);

    // Returns -1 if the name isn't a mol instance name.
    int
    getMolNdxByName(const std::string& rName) const;

    int
    mustGetMolNdxByName(xmlpp::Node* pRequestingNode,
			const std::string& rInstanceName) const
      throw(plx::unknownMolInstanceXcpt);

    // Returns 0 if the name isn't a mol instance name.
    bnd::mol*
    getMolByName(const std::string& rName) const;

    bnd::mol*
    mustGetMolByName(xmlpp::Node* pRequstingNode,
		     const std::string& rInstanceName) const
      throw(plx::unknownMolInstanceXcpt);
  };
}

#endif // PARSERPLEX_H
