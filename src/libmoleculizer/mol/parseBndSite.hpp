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

#ifndef MOL_PARSEBNDSITE_H
#define MOL_PARSEBNDSITE_H

#include "utl/dom.hpp"
#include "mol/mzrBndSite.hpp"
#include "mol/parseSiteShapeName.hpp"
#include "mol/molEltName.hpp"

namespace bnd
{
    class parseMzrBndSite :
        public std::unary_function<xmlpp::Node*, mzrBndSite>
    {
    public:
        mzrBndSite
        operator()( xmlpp::Node* pBindingSiteNode ) const
            throw( utl::xcpt );

    };
}

#endif // MOL_PARSEBNDSITE_H
