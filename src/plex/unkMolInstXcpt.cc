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

#include <sstream>
#include "plex/unkMolInstXcpt.hh"

namespace plx
{
std::string
unkMolInstXcpt::
mkMsg( const std::string& rMolInstName,
       const xmlpp::Node* pOffendingNode )
{
    std::ostringstream msgStream;
    msgStream << utl::dom::xcpt::mkMsg( pOffendingNode )
    << "Unknown mol instance `"
    << rMolInstName
    << ".";
    return msgStream.str();
}

unkMolInstXcpt::
unkMolInstXcpt( const std::string& rMolInstName,
                const xmlpp::Node* pOffendingNode ) :
        utl::xcpt( mkMsg( rMolInstName,
                          pOffendingNode ) )
{}
}
