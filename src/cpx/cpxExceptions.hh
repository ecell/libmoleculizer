/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008 The Molecular Sciences Institute
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Original Author:
//   Nathan Addy, Research Assistant     Voice: 510-981-8748
//   The Molecular Sciences Institute    Email: addy@molsci.org  
//                     
//   
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
