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
#include <boost/foreach.hpp>
#include <vector>
#include <iostream>

#include "../permutation.hh"

using namespace boost::unit_test;
using namespace nmr;
using namespace std;

#define declare_test_suite( test_title ) \
    framework::master_test_suite().p_name.value = test_title;

#define add_test( test_name ) \
    framework::master_test_suite().add( BOOST_TEST_CASE(&test_name) );

#define define_error( type_of_error ) const bool type_of_error( false );

define_error ( EXPECTED_EXCEPTION );
define_error ( DID_NOT_EXPECT_EXCEPTION );

std::vector<int> range (unsigned int max)
{
    std::vector<int> theRange;
    for (unsigned int i = 0;
            i != max;
            ++i)
    {
        theRange.push_back (i);
    }
    return theRange;
}

void permutation_construction_tests()
{
    Permutation nullPerm;
    BOOST_CHECK (nullPerm.getDimension() == 0); // non-critical test => continue after failure

    Permutation perm (10);

    BOOST_REQUIRE (perm.getDimension() == 10 );

    for (int ndx = 0; ndx != perm.getDimension(); ++ndx)
    {
        BOOST_CHECK (perm.getValueAtPosition (ndx) == Permutation::UNDEF );
    }

}


void permutation_value_setting_and_inverse_tests()
{
    Permutation perm (5);
    BOOST_FOREACH (int i, range (5) )
    {
        if (i != 2)
        {
            perm.setValueAtPosition ( i, (3 * i + 27) % 5 );
        }
    }

    BOOST_CHECK ( perm.getValueAtPosition (0) == 2);
    BOOST_CHECK ( perm.getValueAtPosition (1) == 0);
    BOOST_CHECK ( perm.getValueAtPosition (2) == Permutation::UNDEF);
    BOOST_CHECK ( perm.getValueAtPosition (3) == 1);
    BOOST_CHECK ( perm.getValueAtPosition (4) == 4);

    Permutation inversePermutation = perm.invertPermutation();

    BOOST_CHECK ( inversePermutation.getValueAtPosition (0) == 1);
    BOOST_CHECK ( inversePermutation.getValueAtPosition (1) == 3);
    BOOST_CHECK ( inversePermutation.getValueAtPosition (2) == 0);
    BOOST_CHECK ( inversePermutation.getValueAtPosition (3) == Permutation::UNDEF);
    BOOST_CHECK ( inversePermutation.getValueAtPosition (4) == 4);
}


void permutation_composition_tests()
{

    Permutation perm (5);
    BOOST_FOREACH (int i, range (5) )
    {
        if (i != 4)
        {
            perm.setValueAtPosition ( i, (2*i) % 5 );
        }
    }

    Permutation second (5);
    Permutation third (6);

    BOOST_FOREACH (int i, range (5) )
    {
        if (i != 2)
        {
            third.setValueAtPosition ( i, (3 * i + 27) % 5 );
        }
    }


    try
    {
        perm.of (third);
        BOOST_CHECK ( EXPECTED_EXCEPTION  );
    }
    catch (...)
        {}

    try
    {
        third.of ( perm );
        BOOST_CHECK ( EXPECTED_EXCEPTION );
    }
    catch (...)
        {}

    try
    {
        second.of (perm);
        perm.of (second);
    }
    catch (...)
    {
        BOOST_CHECK ( DID_NOT_EXPECT_EXCEPTION );
    }


    {
        Permutation f (3);
        Permutation g (3);

        f.setValueAtPosition (0, 1);
        f.setValueAtPosition (1, 2);
        f.setValueAtPosition (2, 0);

        g.setValueAtPosition (0, 2);
        g.setValueAtPosition (1, 1);
        g.setValueAtPosition (2, 0);

        Permutation h = g.of (f);

        BOOST_CHECK ( h.getValueAtPosition (0) == 1 );
        BOOST_CHECK ( h.getValueAtPosition (1) == 0 );
        BOOST_CHECK ( h.getValueAtPosition (2) == 2 );

        Permutation i = f.of (g);

        BOOST_CHECK ( i.getValueAtPosition (0) == 0);
        BOOST_CHECK ( i.getValueAtPosition (1) == 2);
        BOOST_CHECK ( i.getValueAtPosition (2) == 1);
    }

    {
        Permutation f (3);
        Permutation g (3);

        f.setValueAtPosition (0, 1);
        f.setValueAtPosition (2, 0);

        g.setValueAtPosition (0, 2);
        g.setValueAtPosition (2, 0);

        Permutation h = g.of (f);
        BOOST_CHECK ( h.getValueAtPosition (0) == Permutation::UNDEF );
        BOOST_CHECK ( h.getValueAtPosition (1) == Permutation::UNDEF );
        BOOST_CHECK ( h.getValueAtPosition (2) == 2 );

        Permutation i = f.of (g);
        BOOST_CHECK ( i.getValueAtPosition (0) == 0);
        BOOST_CHECK ( i.getValueAtPosition (1) == Permutation::UNDEF);
        BOOST_CHECK ( i.getValueAtPosition (2) == 1);
    }

}

void check_set_unset_properties()
{
    Permutation f (3);
    BOOST_CHECK ( f.getValueAtPosition (0) == Permutation::UNDEF );

    f.setValueAtPosition (0,0);
    f.setValueAtPosition (1,1);
    BOOST_CHECK ( f.getValueAtPosition (0) == 0 );
    BOOST_CHECK ( f.getValueAtPosition (1) == 1 );

    f.resetValueAtPosition (0);
    BOOST_CHECK ( f.getValueAtPosition (0) == Permutation::UNDEF );
    BOOST_CHECK ( f.getValueAtPosition (1) == 1 );

}

void permutations_must_be_1to1_check()
{
    Permutation f (3);
    f.setValueAtPosition (0,0);

    try
    {
        f.setValueAtPosition (1,0);
        BOOST_CHECK ( EXPECTED_EXCEPTION );
    }
    catch (...)
        {}

}

void check_generalProperties()
{
    Permutation f (3);
    f.setValueAtPosition (0, 0);
    f.setValueAtPosition (1, 1);
    f.setValueAtPosition (2, 2);

    BOOST_CHECK ( f.getIsComplete() );
    BOOST_CHECK ( !f.getIsIncomplete() );

    Permutation g (3);
    g.setValueAtPosition (0, 0);
    g.setValueAtPosition (2, 2);
    BOOST_CHECK ( !g.getIsComplete() );
    BOOST_CHECK ( g.getIsIncomplete() );

}

void check_get_least_value_not_in_permutation()
{
    Permutation f (3);
    BOOST_CHECK ( f.getLeastValueNotInPermutation() == 0);
    f.setValueAtPosition (2, 1);

    BOOST_CHECK ( f.getLeastValueNotInPermutation() == 0);
    f.setValueAtPosition (1, 0);
    BOOST_CHECK ( f.getLeastValueNotInPermutation() == 2);
    f.resetValueAtPosition (1);
    BOOST_CHECK ( f.getLeastValueNotInPermutation() == 0);

}

void check_preimage_unfixeddomain()
{
}

int factorial (int i)
{
    int fact = 1;
    for (int ndx = 1; ndx <= i; ++ndx)
    {
        fact *= ndx;
    }
    return fact;
}

void check_generate_SN()
{
    {
        Permutation::SetOfPermutations aSet;

        cout << "Sizes of S_n..." << endl;
        BOOST_FOREACH ( int i, range (7) )
        {

            if (i == 0) continue;

            aSet.clear();

            Permutation::generate_Sn (aSet, i);

            int size = aSet.size();
            int o_size = factorial (i);

            cout << i << ":\t" << size << "(actual), " << o_size << "(pred)"<< endl;
            BOOST_CHECK ( size == o_size);
        }

        aSet.clear();
        Permutation::generate_Sn (aSet, 3);
        cout << "S_3 = " << endl;
        int i = 0;
        BOOST_FOREACH ( PermutationCref perm, aSet)
        {

            cout << ++i << ":\t" << perm << endl;
        }

    }


    {
        Permutation::SetOfPermutations aSet;
        std::vector<unsigned int> sigVector;

        // sigVector.push_back(1);
        sigVector.push_back (2);
        sigVector.push_back (3);
        sigVector.push_back (1);

        // [0 1 2 3 4 5 6]
        // [0 1 1 2 2 2 3]

        Permutation::generateAllPermutationsMatchingSignature ( aSet, sigVector);
        BOOST_CHECK ( aSet.size() == factorial (2) * factorial (3) );
        BOOST_FOREACH ( const Permutation& perm, aSet)
        {
            cout << perm << endl;
        }
    }
}

void check_factorial()
{
    BOOST_CHECK ( factorial (1) == 1 );
    BOOST_CHECK ( factorial (2) == 2 );
    BOOST_CHECK ( factorial (3) == 6 );
    BOOST_CHECK ( factorial (4) == 24);
}



test_suite*
init_unit_test_suite ( int, char* [] )
{
    declare_test_suite ("Nmr::Permutation Test Suite");

    // register the test case in test tree and specify number of expected failures so
    // this example will pass at runtime. We expect 2 errors: one from failed check and
    // one from memory acces violation
    add_test (check_factorial);
    add_test (permutation_construction_tests);
    add_test (permutation_value_setting_and_inverse_tests);
    add_test (permutation_composition_tests);
    add_test (permutations_must_be_1to1_check);
    add_test (check_set_unset_properties);
    add_test (check_generalProperties);
    add_test (check_get_least_value_not_in_permutation);
    add_test (check_preimage_unfixeddomain);
    add_test (check_generate_SN);

    return 0;
}

// int main()
// {
//     check_generate_SN();
//     return 0;
// }
