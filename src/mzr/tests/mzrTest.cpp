#include <boost/test/included/unit_test.hpp>
using namespace boost::unit_test;


void test_function()
{
    BOOST_CHECK( 0 == 0); // non-critical test => continue after failure
}


//____________________________________________________________________________//


test_suite*
init_unit_test_suite( int, char* [] ) {
    framework::master_test_suite().p_name.value = "MzrUnit_Test";

    // register the test case in test tree and specify number of expected failures so
    // this example will pass at runtime. We expect 2 errors: one from failed check and 
    // one from memory acces violation
    framework::master_test_suite().add( BOOST_TEST_CASE( &test_function ) );

    return 0;
}

//____________________________________________________________________________//

// EOF
