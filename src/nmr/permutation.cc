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

#include <list>

#include <exception>
#include <functional>
#include <iterator>
#include <algorithm>
#include <map>
#include <boost/foreach.hpp>

namespace nmr
{
    const int Permutation::UNDEF = -1;

    Permutation::Permutation(Dimension n) 
        : 
        thePermutation(n, Permutation::UNDEF), //n values of UNDEF
        theDimension( n )
    {
//         // Assert the dimension to be >= 1
//         if (n == 0)
//         {
//             throw GeneralNmrXcpt("Error in constructing Permutation. Dimension cannot have value '0'");
//         }

        // There is only one permutation in S_{1}, the identity function.
        if(n == 1)
        {
            setValueAtPosition(0, 0);
        }
    }

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

    Permutation::Permutation( PermutationCref aPermutation, 
                              PermutationCref bPermutation)
    {
        std::vector<int> corePerm ( aPermutation.getCorePermutation().begin(),
                                    aPermutation.getCorePermutation().end());

        BOOST_FOREACH( int i, bPermutation.getCorePermutation() )
        {
            corePerm.push_back( i + aPermutation.getDimension() );
        }
        
        thePermutation = corePerm;
        theDimension = thePermutation.size();
    }

    Permutation::Permutation(PermutationCref aPermutation, BindingNdx pos, int value) 
        throw(nmr::BadPermutationConstructorXcpt)
        : 
        thePermutation( aPermutation.getCorePermutation() ),
        theDimension( aPermutation.getDimension() )
    {
        if ( aPermutation.getValueAtPosition(pos) != Permutation::UNDEF )
        {
            throw nmr::BadPermutationConstructorXcpt( aPermutation[pos], pos );
        }

        setValueAtPosition(pos, value);
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
            throw nmr::BadPermutationIndexXcpt( getDimension(), pos);
        }
    }

    void 
    Permutation::setValueAtPosition(BindingNdx pos, unsigned int val)
        throw(nmr::BadPermutationIndexXcpt)
    {
        try
        {
            this->thePermutation.at(pos)=val;
        }
        catch(std::out_of_range& e)
        {
            throw nmr::BadPermutationIndexXcpt( getDimension(), pos);
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
            throw nmr::BadPermutationIndexXcpt( getDimension(), pos);
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
    Permutation::getDimension() const
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


    unsigned int 
    Permutation::getLeastValueNotInPermutation() const
        throw(GeneralNmrXcpt)
    {
        if ( getIsComplete() ) throw GeneralNmrXcpt("Logical Error in Permutation::getLeastValueNotInPermutation.  Function was called on a completed permutation.");

        // This function returns the least positive number that is not yet "fixed" 
        // by this partial permutation.
        // For instance, if the permutation is on [0,1,2,3]
        // and 0->2, 1->0 are the only points fixed so far in this partial function 
        // the value to be returned is 1.
       
        // 1. Copy all non-negative values in the Permutation to a new vector.
        UnsignedIntegerVector positiveValues;
        utl::copy_if(thePermutation.begin(),
                     thePermutation.end(),
                     back_inserter(positiveValues),
                     std::bind2nd(std::greater_equal<int>(), 0));

        // 2. Sort the new vector
        std::sort(positiveValues.begin(), 
                  positiveValues.end());


        //3. Iterate through the new vector, returning the first position such that 
        //   value!=position
        for(unsigned int index = 0; 
            index != positiveValues.size();
            ++index)
        {
            if (positiveValues[index] != index) return index;
        }

        //TODO/4 Write a real explanation here.
        //4. I'm too brain-dead to explain this.  It's what the doctor ordered here.  
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
    

    void 
    Permutation::getPreimage(const int& rangeElement,
                             std::set<unsigned int>& refDomainElements) const
    {
        for(unsigned int domainElement = 0;
            domainElement != getDimension();
            ++domainElement)
        {
            if( getValueAtPosition(domainElement) == rangeElement) 
            {
                refDomainElements.insert( domainElement );
            }
        }

        return;
    }

    void 
    Permutation::maximallyExtend()
    {
        if( getNumberOfUndefElements() == 1)
        {

            std::set<unsigned int> domainSubset;
            getUnfixedDomainElements( domainSubset);
            unsigned int singlyRemainingValue = *domainSubset.begin();
            setValueAtPosition(singlyRemainingValue, 
                               getLeastValueNotInPermutation() );
        }
    }

    unsigned int
    Permutation::getNumberOfFixedElements() const
    {
        return getDimension() - getNumberOfUndefElements();
    }

    unsigned int
    Permutation::getNumberOfUndefElements() const
    {
        std::set<unsigned int> domainSubset;
        getPreimage(Permutation::UNDEF,
                    domainSubset);
        return domainSubset.size();
    }

    void 
    Permutation::getUnfixedDomainElements( std::set<unsigned int>& refDomainElements) const
    {
        getPreimage( Permutation::UNDEF, refDomainElements);
    }

    void
    Permutation::generate_Sn( std::set<Permutation>& setOfPermutations, unsigned int N)
    {
        // Illegal dimension.
        if (N == 0) throw 666;

        setOfPermutations.clear();
        
        if ( N == 1)
        {
            std::vector<int> corePermutation(1, 0);
            try
            {
                setOfPermutations.insert( Permutation( std::vector<int>(1, 0) ));
            }
            catch(...)
            {
                int x = 10;
                x += 10;
            }

            return;
        }
        else
        {
            // recursive step.
            generate_Sn( setOfPermutations, N - 1);

            // This could definitely be more efficient if I make set of permutations a ** or 
            // something.
            std::set< Permutation > tmpSet( setOfPermutations.begin(), setOfPermutations.end() );
            setOfPermutations.clear();


            // For example, an element from S_2 might be [1, 0].
            // To Take this element into N=3, we must add the number 2. 
            const int NEXT_RANGE_ELEMENT = N - 1; 

            for( SetOfPermutations::iterator iter = tmpSet.begin();
                 iter != tmpSet.end();
                 ++iter)
            {
                std::list<int> prevDimPerm( iter->getCorePermutation().begin(), 
                                            iter->getCorePermutation().end() );

                std::list<int>::iterator permIter = prevDimPerm.begin();
                
                do
                {
                    permIter = prevDimPerm.insert( permIter, NEXT_RANGE_ELEMENT );
                    setOfPermutations.insert( Permutation( std::vector<int>( prevDimPerm.begin(),
                                                                             prevDimPerm.end())));
                    prevDimPerm.erase( permIter);
                }

                while( permIter++ != prevDimPerm.end());
            }
        }
    }
    


    void
    Permutation::generateAllPermutationsMatchingSignature( SetOfPermutationsRef permSet,
                                                           const std::vector<unsigned int>& signature)
    {
        DECLARE_TYPE( std::vector<unsigned int>, Signature);

        // Make permSet = to S_0
        permSet.clear();
        permSet.insert( Permutation(0) );
        
        for(Signature::const_iterator iter = signature.begin();
            iter != signature.end();
            ++iter)
        {

            // Copy from prev to here.
            SetOfPermutations originalSet(permSet.begin(),
                                          permSet.end());
            permSet.clear();

            // Create the next tmpSet
            SetOfPermutations tmpSet;
            generate_Sn(tmpSet, *iter);

            for( SetOfPermutations::iterator origSetIter = originalSet.begin();
                 origSetIter != originalSet.end();
                 ++origSetIter)
            {

                for( SetOfPermutations::iterator tmpSetIter = tmpSet.begin();
                     tmpSetIter != tmpSet.end();
                     ++tmpSetIter)
                {
                    permSet.insert( Permutation( *origSetIter, *tmpSetIter) );

                }

            }

        }

    }

    void
    Permutation::print() const
    {
        std::cout << "[";
        for(unsigned int ii = 0; ii != thePermutation.size() - 1; ++ii)
        {
            std::cout << thePermutation[ii] << ", ";
        }
        
        std::cout << *thePermutation.rbegin() << ']' << std::endl;

    }

}


