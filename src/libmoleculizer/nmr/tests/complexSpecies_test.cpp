//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2009 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Nathan Addy, Scientific Programmer, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include "nmr/namedMolecule.hpp"
#include "nmr/complexSpecies.hpp"
using namespace boost::unit_test;
using namespace nmr;

#define declare_test_suite( test_title ) \
    framework::master_test_suite().p_name.value = test_title;

#define add_test( test_name ) \
    framework::master_test_suite().add( BOOST_TEST_CASE(&test_name) );

#define define_error( type_of_error ) const bool type_of_error( false );

define_error( EXPECTED_EXCEPTION );
define_error( DID_NOT_EXPECT_EXCEPTION );

std::vector<int> range( unsigned int max )
{
    std::vector<int> theRange;
    for ( unsigned int i = 0;
            i != max;
            ++i )
    {
        theRange.push_back( i );
    }
    return theRange;
}

MolSharedPtr makeMolType_A()
{
    MinimalMolSharedPtr ptrMol( new MinimalMol( "A-Type" ) );
    ptrMol->addNewBindingSite( "One" );
    ptrMol->addNewBindingSite( "Two" );
    ptrMol->addNewBindingSite( "Three" );

    ptrMol->addNewModificationSite( "Site-One", "NULL" );
    return ptrMol;
}

MolSharedPtr makeMolType_B()
{
    MinimalMolSharedPtr ptrMol( new MinimalMol( "B-Type" ) );
    ptrMol->addNewBindingSite( "One" );

    ptrMol->addNewModificationSite( "Site-One", "NULL" );
    return ptrMol;
}


void test_constructors()
{
    ComplexSpecies s;

    BOOST_CHECK( s.getNumberOfMolsInComplex() == 0 );

    MolSharedPtr pMolOne = makeMolType_A();
    MolSharedPtr pMolTwo = makeMolType_A();
    MolSharedPtr pMolThree = makeMolType_B();

    s.addMolToComplex( pMolOne, "Mol_1" );
    s.addMolToComplex( pMolTwo, "Mol_2" );
    s.addMolToComplex( pMolThree, "Mol_3" );

    BOOST_CHECK( s.getNumberOfMolsInComplex() == 3 );
    BOOST_CHECK( s.getNumberOfBindingsInComplex() == 0 );

    s.addBindingToComplex( "Mol_1", "One",
                           "Mol_2", "Two" );

    s.addBindingToComplex( "Mol_2", "One",
                           "Mol_3", "One" );

    BOOST_CHECK( s.getNumberOfBindingsInComplex() == 2 );
}

void test_repr()
{
    ComplexSpecies s;

    MolSharedPtr pMolOne = makeMolType_A();
    MolSharedPtr pMolTwo = makeMolType_A();
    MolSharedPtr pMolThree = makeMolType_B();

    pMolThree->updateModificationState( "Site-One", "NONNULL" );

    s.addMolToComplex( pMolOne, "Mol_1" );
    s.addMolToComplex( pMolTwo, "Mol_2" );
    s.addMolToComplex( pMolThree, "Mol_3" );

    s.addBindingToComplex( "Mol_1", "One",
                           "Mol_2", "Two" );

    s.addBindingToComplex( "Mol_2", "One",
                           "Mol_3", "One" );

    ComplexOutputState aCOS;
    s.constructOutputState( aCOS );

    std::cout << aCOS.repr() << std::endl;

}

test_suite*
init_unit_test_suite( int, char* [] )
{
    declare_test_suite( "Nmr::ComplexSpecies Test Suite" );
    add_test( test_constructors );
    add_test( test_repr );

    return 0;
}

