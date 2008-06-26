#include <boost/test/included/unit_test.hpp>
#include <boost/foreach.hpp>
#include "mzr/moleculizer.hh"
using namespace boost::unit_test;
using namespace mzr;

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

void test_scaffold()
{
    moleculizer theMoleculizer;
    
    // This is pretty bad, as it relies on the program being run from /tests/.libs/.  I have very little 
    // sense of whether this is portable or what.
    theMoleculizer.attachFileName("/home/naddy/Sources/libmoleculizer/src/mzr/tests/scaffold.xml");
}

test_suite*
init_unit_test_suite( int, char* [] ) {
    declare_test_suite( "Moleculizer Test Suite");
    add_test( test_scaffold);

    return 0;
}

