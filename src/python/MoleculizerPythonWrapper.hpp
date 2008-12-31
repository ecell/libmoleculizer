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
//   Nathan Addy, Scientific Programmer, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

#include <vector>
#include <boost/python.hpp>
#include <boost/foreach.hpp>
#include "ReactionNetworkGenerator.hpp"

using namespace boost::python;

template <typename T>
class Convert_Vector_to_Python
{
public:
    static PyObject*
    convert( const std::vector<T>& theVector )
    {
        list myList;
        BOOST_FOREACH( const T& TObject, theVector )
        {
            object myObject( TObject );
            incref( myObject.ptr() );
            myList.append( myObject );
        }
        
        return incref( myList.ptr() );
    }
    
};

#define DECLARE_DEFAULT_VECT_OF_TYPE_TO_PYTHON_CONVERTER( mytype )      \
    to_python_converter< std::vector<mytype>, Convert_Vector_to_Python<mytype> >();
