//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

#ifndef FND_QUERYSPECIESDUMPABLE_H
#define FND_QUERYSPECIESDUMPABLE_H

#include "fnd/query.hh"
#include "fnd/multiSpeciesDumpable.hh"

namespace fnd
{
// This doeesn't include any output (state dump) functionality.  Decided to
// keep xml parsing and output out of template components like this.

template<class speciesT,
class dumpArgT>
class querySpeciesDumpable :
            public multiSpeciesDumpable<speciesT,
            dumpArgT>
{
protected:
    query<speciesT>& rQuery;

public:
// Note that the query is not copied; it's assumed to remain where it is.
    querySpeciesDumpable( const std::string& rName,
                          query<speciesT>& rSpeciesQuery ) :
            multiSpeciesDumpable<speciesT,
            dumpArgT> ( rName ),
            rQuery( rSpeciesQuery )
    {}

    ~querySpeciesDumpable( void )
    {}

// Overrides multiSpeciesDumpable::respond so that only species
// that pass the test get onto the dump list.
    void
    respond( const newSpeciesStimulus<speciesT>& rStimulus )
    {
        const speciesT* pNewSpecies
        = rStimulus.getSpecies();

        if ( rQuery( *pNewSpecies ) )
        {
            this->dumpedSpecies.push_back( pNewSpecies );
        }
    }
};
}

#endif // FND_QUERYSPECIESDUMPABLE_H
