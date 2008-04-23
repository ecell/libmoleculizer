/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2007, 2008  Nathan Addy
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Contact information:
//   Nathan Addy, Research Associate     Voice: 510-981-8748
//   The Molecular Sciences Institute    Email: addy@molsci.org  
//   2168 Shattuck Ave.                  
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#ifndef __NMR_PERMUTATION_HH
#define __NMR_PERMUTATION_HH

#include "nmr/nmrExceptions.hh"

#include "utl/macros.hh"
#include <ostream>
#include <vector>

namespace nmr
{
        // The Permutation class represents a "partial permutation" on Z_n = {0, ..., n-1}.
        // A "partial permutation" is a relation on Z_n X Z_n that can be extended to a 
        // permutation (which is a bijection from Z_n -> Z_n).

        // In this class, this is represented as a function from Z_n -> Z_n U Permutation::UNDEF, 
        // where any domain element that conceptually has no partner in the range Z_n yet is 
        // mapped to Permutation::UNDEF.  The class member function checkPermutationLegality
        // ensures that the Permutation is pre-bijective -- ie that any element of the range which 
        // is not equal to Permutation::UNDEF only shows up once.

        DECLARE_CLASS( Permutation );
        struct Permutation
        {
            // NOTES, TODO:
            // Right now the situation is that the domain is {0, ... , n-1}
            // and the range is {0, ... , n-1} U Permutation::UNDEF.
            // Consequently any of the domain elements (the indexes) are 
            // represented as "unsigned ints" and all the range elements are
            // represented as "ints".  If I ever get the chance, I'd like to come
            // up with better types here.  (although it works along perfectly fine
            // at the moment)

        public:
            static const int UNDEF;

        public:
            DECLARE_TYPE( unsigned int, BindingNdx);
            DECLARE_TYPE( unsigned int, Dimension);
            DECLARE_TYPE( std::vector<int>, IntegerVector);
            DECLARE_TYPE( IntegerVector, CorePermutationType);
            
        public:

            // Constructors 
            //

            // This makes a completely undefined permutation on n objects.
            Permutation(Dimension dim = 0);
            Permutation(PermutationCref aPermutation);
            Permutation(CorePermutationTypeCref aPermutationVector);

            // This constructor copies aPermutation, and then adds the point perm(pos)=value to it.
            Permutation(PermutationCref aPermutation, BindingNdx pos, int value) 
                throw(nmr::BadPermutationConstructorXcpt);
  

            // API
            // 
            int getValueAtPosition(BindingNdx pos) const 
                throw( nmr::BadPermutationIndexXcpt );

            void setValueAtPosition(BindingNdx pos, int val) 
                throw( nmr::BadPermutationIndexXcpt) ;

            void resetValueAtPosition(BindingNdx pos) 
                throw(nmr::BadPermutationIndexXcpt);

            // Returns the Permutation (*this)( compositionPermutation(x) )
            Permutation 
            of(PermutationCref compositionPermutation) const
                throw( nmr::IncompatiblePermutationsXcpt);

            Permutation 
            invertPermutation() const;
  
            bool
            getIsBijection() const;

            bool 
            getIsComplete() const;

            bool 
            getIsIncomplete() const ;

            Dimension 
            getPermutationSize() const;
            
            Dimension
            getDimension() const
            {
                return getPermutationSize();
            }
            
            int& operator[](const BindingNdx& n) 
                throw(nmr::BadPermutationIndexXcpt);

            const int& operator[](const BindingNdx& n) const
                throw(nmr::BadPermutationIndexXcpt);

            unsigned int 
            getLeastValueNotInPermutation() const 
                throw(nmr::GeneralNmrXcpt);

            bool 
            checkPermutationLegality() const;

            bool 
            operator==(PermutationCref pm);

            bool 
            operator<(PermutationCref pm) const;

            CorePermutationTypeCref
            getCorePermutation() const
            {
                return thePermutation;
            }

        protected:
            CorePermutationType thePermutation;
            Dimension theDimension;
        };

}

#endif

