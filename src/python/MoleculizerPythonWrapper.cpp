#include <boost/python.hpp>
#include "mzr/moleculizer.hh"
using namespace boost::python;
using namespace mzr;

BOOST_PYTHON_MODULE( Moleculizer )
{
  class_<moleculizer>("ReactionNetworkGenerator")
      .def("greeting", &moleculizer::greeting);
    
}
