#include <boost/python.hp>
using namespace boost::python;
using namespace mzr;

BOOST_PYTHON_MODULE( Moleculizer )
{
  class_<moleculizer>("ReactionNetworkGenerator")
    .def(attachFile, &moleculizer::attachFileName)
    .def(attachString, &moleculizer:attachString)
    .def(findUnaryReaction, &moleculizer::findReactionWithSubstrates
    
}
