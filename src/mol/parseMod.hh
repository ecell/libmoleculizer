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

#ifndef MOL_PARSEMOD_H
#define MOL_PARSEMOD_H

#include "utl/dom.hh"
#include "cpx/modification.hh"
#include "mol/molEltName.hh"

namespace bnd
{
class parseModification :
            public std::unary_function<xmlpp::Node*, cpx::modification*>
{
public:
    cpx::modification*
    operator()( xmlpp::Node* pModificationNode ) const
    {
        xmlpp::Element* pModElt
        = utl::dom::mustBeElementPtr( pModificationNode );

        std::string name
        = utl::dom::mustGetAttrString( pModElt,
                                       eltName::modification_nameAttr );

        xmlpp::Element* pWeightDeltaElt
        = utl::dom::mustGetUniqueChild( pModElt,
                                        eltName::weightDelta );

        double weightDelta = utl::dom::mustGetAttrDouble
                             ( pWeightDeltaElt,
                               eltName::weightDelta_daltonsAttr );

        return new cpx::modification( name,
                                      weightDelta );
    }
};
}

#endif // MOL_PARSEMOD_H
