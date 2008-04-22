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

#include "nmr/permutation.hh"
#include "utl/utility.hh"

#include <exception>
#include <functional>
#include <iterator>
#include <algorithm>
#include <map>

namespace nmr
{
    const int Permutation::UNDEF = -1;

        Permutation::Permutation(Dimension n) 
            : 
            thePermutation(n, Permutation::UNDEF), //n values of UNDEF
            theDimension( n )
        {}

        Permutation::Permutation(PermutationCref aPermutation) 
            : 
            thePermutation( aPermutation.getCorePermutation() ),
            theDimension( aPermutation.getDimension() )
        {}

        Permutation::Permutation(CorePermutationTypeCref aPermutationVector) 
            : 
            thePermutation(aPermutationVector.begin(), 
                           aPermutationVector.end()),
            theDimension( aPermutationVector.size() )
        {}

        Permutation::Permutation(PermutationCref aPermutation, BindingNdx pos, int value) 
            throw(nmr::BadPermutationConstructorXcpt)
            : 
            thePermutation(aPermutation.getCorePermutation() ),
            theDimension( aPermutation.getPermutationSize() )
        {
            if ( aPermutation.getValueAtPosition(pos) != Permutation::UNDEF )
            {
                throw nmr::BadPermutationConstructorXcpt( aPermutation[pos], pos );
            }

            thePermutation[pos]=value;
        }
  
        int 
        Permutation::getValueAtPosition(BindingNdx pos) const
            throw(nmr::BadPermutationIndexXcpt)
        {
            try
            {
                return this->thePermutation.at(pos);
            }
            catch(std::out_of_range& e)
            {
                throw nmr::BadPermutationIndexXcpt( getPermutationSize(), pos);
            }
        }

        void 
        Permutation::setValueAtPosition(BindingNdx pos, int val)
            throw(nmr::BadPermutationIndexXcpt)
        {
            try
            {
                this->thePermutation.at(pos)=val;
            }
            catch(std::out_of_range& e)
            {
                throw nmr::BadPermutationIndexXcpt( getPermutationSize(), pos);
            }
        }

        void 
        Permutation::resetValueAtPosition(BindingNdx pos) 
            throw(nmr::BadPermutationIndexXcpt)
        {
            try
            {
                this->thePermutation.at(pos)= Permutation::UNDEF;
            }
            catch(std::out_of_range& e)
            {
                throw nmr::BadPermutationIndexXcpt( getPermutationSize(), pos);
            }
        }

        Permutation 
        Permutation::of(PermutationCref compositionPermutation) const
            throw( nmr::IncompatiblePermutationsXcpt)
        {
            // If the two dimensions don't match up, someone fucked up.  
            // Let's judge that blockhead by taking exception!
            if( getDimension() != compositionPermutation.getDimension() ) 
            {
                throw nmr::IncompatiblePermutationsXcpt( getDimension(), compositionPermutation.getDimension() );
            }
    
            Permutation tmpPerm( this->getDimension() );
        
            for(BindingNdx i = 0;
                i != this->getDimension();
                ++i)
            {
                int intermediateValue=compositionPermutation.getValueAtPosition(i);
                if( intermediateValue==Permutation::UNDEF )
                {
                    tmpPerm[i]=Permutation::UNDEF;
                }
                else
                {
                    int finalValue=this->getValueAtPosition(intermediateValue);
                    tmpPerm[i]=finalValue;
                }

            }
            
            return tmpPerm;
        }

        Permutation 
        Permutation::invertPermutation() const
        {
            // This function returns a brand new Permutation that is the inverse
            // to (*this).
         
            Permutation invertedPermutation( getDimension() );;

            for( BindingNdx index = 0;
                 index != theDimension;
                 ++index)
            {
                const int value( getValueAtPosition( index ) );
                if ( value != Permutation::UNDEF )
                {
                    invertedPermutation[ value ] = index;
                }
            }

            return invertedPermutation;
        }


        bool 
        Permutation::getIsComplete() const
        {
            // A permutation is complete iff it has no undefinded values and 
            // is also a legal permutation.
            CorePermutationType::const_iterator i=find(thePermutation.begin(),
                                                       thePermutation.end(),
                                                       Permutation::UNDEF);
            return (checkPermutationLegality() && i == thePermutation.end() );
        }

        bool
        Permutation::getIsBijection() const
        {
            return getIsComplete();
        }


        bool 
        Permutation::getIsIncomplete() const
        {
            return !getIsComplete();
        }

        Permutation::Dimension 
        Permutation::getPermutationSize() const
        {
            return theDimension;
        }

        int& 
        Permutation::operator[](const BindingNdx& n)
            throw(nmr::BadPermutationIndexXcpt)
        {
            if (n >= getDimension() ) throw nmr::BadPermutationIndexXcpt( getDimension(), n);

            return thePermutation[n];  
        }

        const int& 
        Permutation::operator[](const BindingNdx& n) const
            throw(nmr::BadPermutationIndexXcpt)
        {
            if (n >= getDimension() ) throw nmr::BadPermutationIndexXcpt( getDimension(), n);

            return thePermutation[n];  
        }


        int 
        Permutation::getLeastValueNotInPermutation() const
        {
            //TODO: Check the accuracy of this function.

            //copy the Permutation to a new vector, ommitting any element where the value is less than 0
            //sort the new vector
            //iterate through the new vector, returning the first position such that value!=position
            IntegerVector positiveValues;

            utl::copy_if(thePermutation.begin(),
                         thePermutation.end(),
                         back_inserter(positiveValues),
                         std::bind2nd(std::greater_equal<int>(), 0));

            std::sort(positiveValues.begin(), 
                      positiveValues.end());

            for(int i=0; i!=(int) positiveValues.size();++i)
            {
                if (positiveValues[i] != i)
                {
                    return i;
                }
            }
            return positiveValues.size();
        }


        bool 
        Permutation::checkPermutationLegality() const
        {

            // A permutation is legal iff
            // 1. Every member if the range is in {0, ..., dimension -1 } U Permutation::UNDEF
            // and 
            // 2. If Permutation(x) == z and Permutation(y) == z where z is an element of 
            //    {0, ..., dim - 1} then x == z.

            IntegerVector rangeCounts(this->thePermutation.size(), 0);

            for(CorePermutationType::const_iterator iter = thePermutation.begin();
                iter !=thePermutation.end();
                ++iter )
            {
                //if (*iter) isn't in -1, 0, 1,..., thePermutation.size()-1, return false
                if( (*iter)<Permutation::UNDEF || *iter >= static_cast<int>(getDimension()) )
                {
                    return false;
                }
                else if( *iter != Permutation::UNDEF)
                {
                    rangeCounts[*iter] += 1;
                }
            }

            for(std::vector<int>::const_iterator iter =rangeCounts.begin();
                iter != rangeCounts.end();
                ++iter)
            {
                if ( *iter > 1)
                {
                    return false;
                }
            }

            return true;
        }


        bool 
        Permutation::operator==(const Permutation& pm)
        {
            if ((this->thePermutation)==(pm.thePermutation))
                return true;
            else return false;
        }

        bool 
        Permutation::operator<(const Permutation& pm) const
        {

            // TODO: is this correct?
            if ( this->getDimension() < pm.getDimension() ) return true;
            else if ( this->getDimension() > pm.getDimension() ) return false;

            CorePermutationType::const_iterator jjIter;
            for(CorePermutationType::const_iterator iiIter = thePermutation.begin(), 
                    jjIter = pm.getCorePermutation().begin();
                iiIter != thePermutation.end();
                ++iiIter, ++jjIter)
            {
                if (*iiIter < *jjIter ) return true;
            }

            return false;
        }
}


