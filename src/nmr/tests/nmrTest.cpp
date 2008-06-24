#include <boost/test/included/unit_test.hpp>
#include "../permutation.hh"
using namespace boost::unit_test;


void permutation_test_function()
{
    nmr::Permutation pp;
    
    BOOST_CHECK(pp.getDimension() == 0); // non-critical test => continue after failure

}



//____________________________________________________________________________//


test_suite*
init_unit_test_suite( int, char* [] ) {
    framework::master_test_suite().p_name.value = "FooTest";

    // register the test case in test tree and specify number of expected failures so
    // this example will pass at runtime. We expect 2 errors: one from failed check and 
    // one from memory acces violation
    framework::master_test_suite().add( BOOST_TEST_CASE( &permutation_test_function ) );

    return 0;
}

//____________________________________________________________________________//

// EOF
