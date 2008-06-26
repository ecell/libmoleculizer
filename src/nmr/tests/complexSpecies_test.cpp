#include <boost/test/included/unit_test.hpp>
#include "nmr/namedMolecule.hh"
#include "nmr/complexSpecies.hh"
using namespace boost::unit_test;
using namespace nmr;

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

MolSharedPtr makeMolType_A()
{
    MinimalMolSharedPtr ptrMol( new MinimalMol("A-Type") );
    ptrMol->addNewBindingSite("One");
    ptrMol->addNewBindingSite("Two");
    ptrMol->addNewBindingSite("Three");
    
    ptrMol->addNewModificationSite("Site-One", "NULL");
    return ptrMol;
}

MolSharedPtr makeMolType_B()
{
    MinimalMolSharedPtr ptrMol( new MinimalMol("A-Type") );
    ptrMol->addNewBindingSite("One");

    ptrMol->addNewModificationSite("Site-One", "NULL");
    return ptrMol;
}


void test_constructors()
{
    ComplexSpecies s;

    BOOST_CHECK( s.getNumberOfMolsInComplex() == 0);

    MolSharedPtr pMolOne = makeMolType_A();
    MolSharedPtr pMolTwo = makeMolType_A();
    MolSharedPtr pMolThree = makeMolType_B();

    s.addMolToComplex(pMolOne, "Mol_1");
    s.addMolToComplex(pMolTwo, "Mol_2");
    s.addMolToComplex(pMolThree, "Mol_3");
    
    BOOST_CHECK( s.getNumberOfMolsInComplex() == 3);
    BOOST_CHECK( s.getNumberOfBindingsInComplex() == 0);

    s.addBindingToComplex("Mol_1", "One", 
                          "Mol_2", "Two");

    s.addBindingToComplex("Mol_2", "One", 
                          "Mol_3", "One");

    BOOST_CHECK( s.getNumberOfBindingsInComplex() == 2);
}

test_suite*
init_unit_test_suite( int, char* [] ) {
    declare_test_suite( "Nmr::ComplexSpecies Test Suite");
    add_test(test_constructors);

    return 0;
}

