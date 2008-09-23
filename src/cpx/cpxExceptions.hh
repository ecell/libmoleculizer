/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008 The Molecular Sciences Institute
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
//   Nathan Addy, Research Assistant     Voice: 510-981-8748
//   The Molecular Sciences Institute    Email: addy@molsci.org  
//   2168 Shattuck Ave.                  
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#ifndef CPX_EXCEPTIONS_H
#define CPX_EXCEPTIONS_H

#include "utl/xcpt.hh"

namespace cpx
{

    class plexIsNotSimpleGraphXcpt;



    class plexIsNotSimpleGraphXcpt :
        public utl::xcpt
    {
        static std::string
        mkMsg();

    public:
        plexIsNotSimpleGraphXcpt(void)
            :
            utl::xcpt(mkMsg())
        {}
    };
}

#endif
