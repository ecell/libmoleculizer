//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2009 The Molecular Sciences Institute.
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

#ifndef CPX_ISOSEARCH_H
#define CPX_ISOSEARCH_H

#include "cpx/plexIso.hh"

namespace cpx
{
    template<class plexT>
    class isoSearch
    {
        const plexT& rLeft;
        const plexT& rRight;
        
        // Determines if rCurrentIso can be extended over all the bindings
        // starting at leftBindingIndex in the left plex.  This is the basic
        // recursive step in the process of finding an injection or isomorphism.
        bool
        mapRestBindings( int leftBindingIndex,
                         const plexIso& rCurrentIso ) const;
        
        
    public:
        isoSearch( const plexT& rLeftPlex,
                   const plexT& rRightPlex ) :
            rLeft( rLeftPlex ),
            rRight( rRightPlex )
        {}
        
        virtual
        ~isoSearch( void )
        {}
        
        // Virtual function called on the first isomorphism (or injection) found
        // during the search, if any.  This could be to copy the isomorphism out
        // as a return value, as in the derived reportIsoSearch class.
        virtual void
        onSuccess( const plexIso& rIso ) const
        {}
        
        // Determines if the left plex is a subcomplex of the right plex.
        bool
        findInjection( void ) const;
        
        // Determines if the left plex is isomorphic to the right plex.
        bool
        findIso( void ) const;
    };
}

#include "cpx/isoSearchImpl.hh"

#endif // CPX_ISOSEARCH_H
