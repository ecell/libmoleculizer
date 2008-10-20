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

#ifndef __PROPERTIEDCLASS_HPP
#define __PROPERTIEDCLASS_HPP

#include "fnd/utility.hh"
#include "propXcpt.hpp"
#include "propertyValue.hpp"
#include "property.hpp"


class PropertiedClass
{
public:
    PropertiedClass();
    PropertiedClass( const PropertiedClass& propertiedClass );
    ~PropertiedClass();

    void
    addProperty( Property propertyName );

    void
    createProperty( String propertyName );

    Value operator[]( const String& propertyName ) throw( PropertyDoesNotExistXcpt );
    Value getPropertyValue( const String& propertyName ) throw( PropertyDoesNotExistXcpt );

    Property& getProperty( const String& propertyName );
    const Property getProperty( const String& propertyName ) const;

protected:

    std::map<const String*, Property, fnd::aux::compareByPtrValue<String> > thePropertyMap;
};


#endif
