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


Value::Value()
    :
    underlyingValue( new Real(0.0f) )
{}

Value::~Value()
{}

void Value::setValue(const Value& value)
{
    *underlyingValue = value->valueAsReal();
}

void Value::setValue(String strValue)
{
    throw utl::NotImplementedYet("void Value::setValue(String strValue)");
}

void Value::setValue(Integer intValue)
{
    *underlyingValue = intValue;
}

void Value::setValue(Real realValue)
{
    *underlyingValue = value->valueAsReal();
}

String Value::valueAsString() const
{
    throw utl::NotImplementedYet("String Value::valueAsString()");
}

Integer Value::valueAsInteger() const
{
    return static_cast<Integer>( this->valueAsReal );
}

Real Value::valueAsReal() const
{
    return *underlyingValue;
}

Value Value::clone() const
{
    return Value( this->valueAsReal() );
}
