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
    convert( const std::vector<T>& theVector)
    {
        list myList;
        BOOST_FOREACH(const T& TObject, theVector)
        {
            object myObject( TObject );
            incref( myObject.ptr() );
            myList.append( myObject );
        }
        
        return incref( myList.ptr() );
    }
    
};

#define DECLARE_DEFAULT_VECT_OF_TYPE_TO_PYTHON_CONVERTER( mytype )\
    to_python_converter< std::vector<mytype>, Convert_Vector_to_Python<mytype> >();
