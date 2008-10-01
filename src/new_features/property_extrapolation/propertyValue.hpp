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
//   Nathan Addy, Scientific Programmer, Molecular Sciences Institute
//
// Modifing Authors:
//
//

#ifndef __VALUE_HPP
#define __VALUE_HPP

#include "utl/defs.hpp" 
#include "utl/common.hpp"

using namespace boost;

// This class represents a polymorphic value, something that can be used
// to hold a real, integer, or string.  That is, ultimately it will be that.
// For now, values must be Reals.


    


class Value
{
public:
    Value();
    Value(const Value& value);
    ~Value();

    const Value& operator=(const Value& value);

    // Initially, I am planning to have the different Value types do 
    // template overriding, and so I am not making these virtual.  Should I?
    void setValue(const Value& value);
    void setValue(String strValue) throw( utl::NotImplementedXcpt );
    void setValue(Integer intValue);
    void setValue(Real realValue);

    String valueAsString() const;
    Integer valueAsInteger() const;
    Real valueAsReal() const;
    
    ValuePtr clone() const;

protected:
    typedef shared_ptr<Real> UnderlyingValuePtr;
    UnderlyingValuePtr  underlingValue;
};

#include "propertyValueImpl.hpp"

#endif
