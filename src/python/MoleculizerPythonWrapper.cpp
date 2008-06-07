#include <vector>
#include <boost/python.hpp>
#include "ReactionNetworkGenerator.hpp"
using namespace boost::python;


class ReactionVector_to_Python
{
public:
    static PyObject* convert( const std::vector<Reaction>& theReactionVector );
};


class StringVector_to_Python
{
public:
    static PyObject* convert(const std::vector<std::string>& theStringVector);
};


BOOST_PYTHON_MODULE( Moleculizer )
{
    class_<Reaction>("Reaction")
        .def("getRate", &Reaction::getRate)
        .def("getSubstrates", &Reaction::getSubstrates)
        .def("getProducts", &Reaction::getProducts)
        .add_property("rate", &Reaction::getRate)
        .add_property("substrates", &Reaction::getSubstrates)
        .add_property("products", &Reaction::getProducts)
        ;

    class_<ReactionNetworkGenerator>("ReactionNetworkGenerator")
        .def("DEBUG", &ReactionNetworkGenerator::DEBUG)
        .def("addRules", &ReactionNetworkGenerator::addRules)
        .def("getBinaryReactions", &ReactionNetworkGenerator::getBinaryReactions)
        .def("getUnaryReactions", &ReactionNetworkGenerator::getUnaryReactions)
        .def("showAllSpecies", &ReactionNetworkGenerator::showAllSpecies)
        .def("checkSpeciesExists", &ReactionNetworkGenerator::checkSpeciesExists)
        .def("checkSpeciesNameLegality", &ReactionNetworkGenerator::checkSpeciesNameLegality)
        ;

    to_python_converter< std::vector<Reaction>, ReactionVector_to_Python>();
    to_python_converter< std::vector<std::string>, StringVector_to_Python>();
}


PyObject*
ReactionVector_to_Python::convert(const std::vector<Reaction>& theReactionVector)
{
    list myList;
    
    for(unsigned int i = 0; i!= theReactionVector.size(); ++i)
    {
        object myObject( theReactionVector[i] );

        incref( myObject.ptr() );

        myList.append(myObject);
    }
    
    return incref( myList.ptr() );
}


PyObject*
StringVector_to_Python::convert(const std::vector<std::string>& theStringVector)
{
    list myList;
    for(unsigned int i = 0; i != theStringVector.size(); ++i)
    {
        object myObject( theStringVector[i]);
        incref( myObject.ptr() );
        
        myList.append( myObject);
    }

    return incref( myList.ptr() );
}
