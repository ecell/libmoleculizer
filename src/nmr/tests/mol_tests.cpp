//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                                                                          
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
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
#include "nmr/namedMolecule.hh"

#ifdef HAVE_CONFIG_H
#include "moleculizer_config.hh"
#endif

using namespace boost::unit_test;
using namespace nmr;
using namespace std;

#define declare_test_suite( test_title ) \
    framework::master_test_suite().p_name.value = test_title;

#define add_test( test_name ) \
    framework::master_test_suite().add( BOOST_TEST_CASE(&test_name) );

#define define_error( type_of_error ) const bool type_of_error( false );

define_error( EXPECTED_EXCEPTION_HERE );
define_error( DID_NOT_EXPECT_EXCEPTION_HERE );

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

void test_constructors()
{
    MinimalMol mm("DebugMolType");
    MinimalMol mm_copy( mm );

    BOOST_CHECK( mm.getMolType() == "DebugMolType");
    BOOST_CHECK( mm_copy.getMolType() == "DebugMolType");
}

void test_semantics()
{

    MinimalMol mm("FooType");
    mm.addNewBindingSite("C");
    mm.addNewBindingSite("A");
    mm.addNewBindingSite("B");
    mm.addNewBindingSite("D");
        
    mm.addNewModificationSite("Site-1", "NULL");
    mm.addNewModificationSite("Site-2", "NULL");
    mm.addNewModificationSite("Site-3", "NULL");

    BOOST_CHECK( mm.getBindingSiteInteger("A") == 1);
    BOOST_CHECK( mm.getBindingSiteInteger("C") == 0);


    BOOST_CHECK( mm.checkIfBindingSiteExists("A") );
    BOOST_CHECK( mm.checkIfBindingSiteExists("B") );
    BOOST_CHECK( mm.checkIfBindingSiteExists("C") );
    BOOST_CHECK( mm.checkIfBindingSiteExists("D") );
    
    BOOST_CHECK( mm.checkIfBindingSiteExists("E") == false );

    BOOST_CHECK( mm.checkIfModificationSiteExists("Site-3") );
    BOOST_CHECK( mm.checkIfModificationSiteExists("Site-1") );
    BOOST_CHECK( mm.checkIfModificationSiteExists("Site-2") );


    BOOST_CHECK( mm.checkIfModificationSiteExists("Site-4") == false );


    BOOST_CHECK( mm.getBindingSiteInteger("C") == 0);
    BOOST_CHECK( mm.getBindingSiteInteger("A") == 1);
    BOOST_CHECK( mm.getBindingSiteInteger("D") == 3);
    BOOST_CHECK( mm.getBindingSiteInteger("B") == 2);

    // BOOST_CHECK( mm.getBindingSiteInteger("A") == 1);
//     BOOST_CHECK( mm.getBindingSiteInteger("B") == 2);
//     BOOST_CHECK( mm.getBindingSiteInteger("C") == 3);
//     BOOST_CHECK( mm.getBindingSiteInteger("D") == 4);

  //   try
//     {
//         mm.getBindingSiteInteger("E");
//         BOOST_CHECK( EXPECTED_EXCEPTION_HERE );
//     }
//     catch(...)
//     {}


    BOOST_CHECK( mm.getModificationSiteInteger("Site-1") == 0);
    BOOST_CHECK( mm.getModificationSiteInteger("Site-2") == 1);
    BOOST_CHECK( mm.getModificationSiteInteger("Site-3") == 2);

    BOOST_CHECK( mm.checkIfBindingSiteIsBound("C") == false );

    mm.bindAtBindingSite("C");

    BOOST_CHECK( mm.checkIfBindingSiteIsBound("C") == true );

    mm.unbindAtBindingSite( "C" );

    BOOST_CHECK( mm.checkIfBindingSiteIsBound("C") == false );

    BOOST_CHECK(mm.getModificationValueAtModificationSite("Site-2") == "NULL");
    mm.updateModificationState("Site-1", "ALPHA");
    mm.updateModificationState("Site-3", "BETA");

    BOOST_CHECK(mm.getModificationValueAtModificationSite("Site-1") == "ALPHA");
    BOOST_CHECK(mm.getModificationValueAtModificationSite("Site-2") == "NULL");
    BOOST_CHECK(mm.getModificationValueAtModificationSite("Site-3") == "BETA");

    mm.updateModificationState("Site-1", "NULL");
    
    BOOST_CHECK(mm.getModificationValueAtModificationSite("Site-1") == "NULL");
    BOOST_CHECK(mm.getModificationValueAtModificationSite("Site-2") == "NULL");
    BOOST_CHECK(mm.getModificationValueAtModificationSite("Site-3") == "BETA");

}

test_suite*
init_unit_test_suite( int, char* [] ) {
    declare_test_suite( "nmr::Mol Test Suite");

    add_test( test_constructors );
    add_test( test_semantics );

    return 0;
}

