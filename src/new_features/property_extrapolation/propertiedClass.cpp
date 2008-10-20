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

#include "propertiedClass.hpp"

PropertiedClass::PropertiedClass()
{}


PropertiedClass::PropertiedClass( const PropertiedClass& propertiedClass )
{

}


PropertiedClass::~PropertiedClass()
{
}

Value
PropertiedClass::operator[]( const String& propertiedName ) throw( PropertyDoesNotExistXcpt )
{
    return this->getPropertyValue( propertiedName );
}

Value
PropertiedClass::getPropertyValue( const String& propertyName ) throw( PropertyDoesNotExistXcpt )
{
    return getProperty( propertyName ).getValue();
}

Property&
PropertiedClass::getProperty( const String& propertyName )
{
    return thePropertyMap.find( &propertyName )->second;
}



void PropertiedClass::addProperty( Property propertyToAdd )
{

    if ( propertyToAdd.hasOwner() )
    {
        throw utl::xcpt( "property '" + propertyToAdd.getName() + "' already has an owner." );
    }

    const String* ptrString = &propertyToAdd.getName();

    if ( thePropertyMap.find( ptrString ) == thePropertyMap.end() )
    {
        propertyToAdd.setOwner( this );
        thePropertyMap.insert( std::make_pair( ptrString, propertyToAdd ) );
    }
    else
    {
        String badPropertyName( *ptrString );
        // throw PropertyAlreadyExistsXcpt( propertyName )
        throw utl::xcpt( "Property '" + badPropertyName + "' already exists in this PropertiedClass." );
    }
}

void PropertiedClass::createProperty( String propertyName )
{

    Property newProp( propertyName, this );

    const String* ptrString = &newProp.getName();

    if ( thePropertyMap.find( ptrString ) == thePropertyMap.end() )
    {
        thePropertyMap.insert( std::make_pair( ptrString, newProp ) );
    }
    else
    {
        // throw PropertyAlreadyExistsXcpt( propertyName )
        throw utl::xcpt( "Property '" + propertyName + "' already exists in this PropertiedClass." );
    }

}
