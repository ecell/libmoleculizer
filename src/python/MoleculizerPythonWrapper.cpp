#include "MoleculizerPythonWrapper.hpp"

BOOST_PYTHON_MODULE( Moleculizer )
{
    class_<Species>("Species", init<const mzr::mzrSpecies&>())
        .add_property("name", &Species::getName)
        .add_property("mass", &Species::getMass)
        ;

    class_<Reaction>("Reaction", init<const Reaction::CoreRxnType&>())
        .def("getRate", &Reaction::getRate)
        .def("getSubstrates", &Reaction::getSubstrates)
        .def("getProducts", &Reaction::getProducts)
        .add_property("rate", &Reaction::getRate)
        .add_property("substrates", &Reaction::getSubstrates)
        .add_property("products", &Reaction::getProducts)
        ;

    class_<BasicComplexRepresentation>("BasicComplexRepresentation")
        .def("addMolToComplex", &BasicComplexRepresentation::addMolNameToComplex)
        .def("addModificationToComplex", &BasicComplexRepresentation::addModificationToComplex)
        .def("addBindingToComplex", &BasicComplexRepresentation::addBindingToComplex)
        ;

    class_<ReactionNetworkGenerator>("ReactionNetworkGenerator")
        .def("runDebugMode", &ReactionNetworkGenerator::runInteractiveMode)
        .def("addRules", &ReactionNetworkGenerator::addRules)
        .def("getBinaryReactions", &ReactionNetworkGenerator::getBinaryReactions)
        .def("getUnaryReactions", &ReactionNetworkGenerator::getUnaryReactions)
        .def ("getSpecies", &ReactionNetworkGenerator::getSpecies)
        .def("checkSpeciesNameLegality", &ReactionNetworkGenerator::checkSpeciesNameLegality)
        .def("showAllSpecies", &ReactionNetworkGenerator::showAllSpecies)
        .def("incrementSpecies", &ReactionNetworkGenerator::incrementSpecies)
        .def("getNumberOfSpecies", &ReactionNetworkGenerator::getNumberOfSpecies)
        .def("getNumberOfReactions", &ReactionNetworkGenerator::getNumberOfReactions)
        .def("showLiveSpecies", &ReactionNetworkGenerator::showLiveSpecies)
        .def("showDeadSpecies", &ReactionNetworkGenerator::showDeadSpecies)
        .def("getNumUnaryReactions",&ReactionNetworkGenerator::getNumUnary)
        .def("getNumBinary", &ReactionNetworkGenerator::getNumBinary)
        .def("generateNameFromBasicComplexRepresentation", &ReactionNetworkGenerator::generateNameFromBasicComplexRepresentationStrict)
        .def("printRxns", &ReactionNetworkGenerator::showAllReactions)
        ;


    DECLARE_DEFAULT_VECT_OF_TYPE_TO_PYTHON_CONVERTER( Reaction );
    DECLARE_DEFAULT_VECT_OF_TYPE_TO_PYTHON_CONVERTER( std::string );
    DECLARE_DEFAULT_VECT_OF_TYPE_TO_PYTHON_CONVERTER( Species );
}

