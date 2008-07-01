#include <boost/test/included/unit_test.hpp>
#include <boost/foreach.hpp>
#include <iostream>
using namespace boost::unit_test;
using std::cout;
using std::endl;


#define declare_test_suite( test_title ) \
    framework::master_test_suite().p_name.value = test_title;

#define add_test( test_name ) \
    framework::master_test_suite().add( BOOST_TEST_CASE(&test_name) );

#define define_error( type_of_error ) const bool type_of_error( false );

define_error( EXPECTED_EXCEPTION );
define_error( DID_NOT_EXPECT_EXCEPTION );

std::vector<int> range(unsigned int max)
{
    std::vector<int> theRange;
    for(unsigned int i = 0;
        i != max;
        ++i)
    {
        theRange.push_back(i);
    }
    return theRange;
}

void test_function()
{
    throw 0;
    BOOST_CHECK( 0 == 0); // non-critical test => continue after failure
}

test_suite*
init_unit_test_suite( int, char* [] ) {
    declare_test_suite( "Foo Test Suite");
    add_test( test_function);

    return 0;
}

