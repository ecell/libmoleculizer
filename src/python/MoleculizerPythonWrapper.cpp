#include "MoleculizerPythonWrapper.hpp"

BOOST_PYTHON_MODULE( Moleculizer )
{
    class_<Species>("Species")
        .add_property("name", &Species::getName)
        .add_property("mass", &Species::getMass)
        ;

    class_<Reaction>("Reaction")
        .def("getRate", &Reaction::getRate)
        .def("getSubstrates", &Reaction::getSubstrates)
        .def("getProducts", &Reaction::getProducts)
        .add_property("rate", &Reaction::getRate)
        .add_property("substrates", &Reaction::getSubstrates)
        .add_property("products", &Reaction::getProducts)
        ;

    class_<ReactionNetworkGenerator>("ReactionNetworkGenerator")
        .def("addRules", &ReactionNetworkGenerator::addRules)
        .def("getBinaryReactions", &ReactionNetworkGenerator::getBinaryReactions)
        .def("getUnaryReactions", &ReactionNetworkGenerator::getUnaryReactions)
        .def ("getSpecies", &ReactionNetworkGenerator::getSpecies)
        .def("checkSpeciesNameLegality", &ReactionNetworkGenerator::checkSpeciesNameLegality)
        .def("showAllSpecies", &ReactionNetworkGenerator::showAllSpecies)
        ;

    DECLARE_DEFAULT_VECT_OF_TYPE_TO_PYTHON_CONVERTER( Reaction );
    DECLARE_DEFAULT_VECT_OF_TYPE_TO_PYTHON_CONVERTER( std::string );
    DECLARE_DEFAULT_VECT_OF_TYPE_TO_PYTHON_CONVERTER( Species );
}

