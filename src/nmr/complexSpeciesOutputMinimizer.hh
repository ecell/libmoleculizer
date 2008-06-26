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

#ifndef __CANONICALNMR_HH
#define __CANONICALNMR_HH

#include "complexSpecies.hh"
#include "permutation.hh"
#include "partialTokenList.hh"
#include "permutationName.hh"
#include "abstractMol.hh"

#include <set>
#include <map>
#include <functional>

namespace nmr
{
    DECLARE_CLASS( BindingSorterByNdx );
    DECLARE_CLASS(ComplexSpeciesOutputMinimizer);
    class ComplexSpeciesOutputMinimizer
    {

    public:

        ComplexOutputState 
        getMinimalOutputState(ComplexSpecies aComplexSpecies);

        ComplexOutputState
        getMinimalOutputStateViaBruteForce( ComplexSpecies aComplexSpecies);

    protected:
        typedef std::pair<int, int> __Range;
        typedef std::pair<int, int> __BindingSite;
        typedef std::pair<__BindingSite, __BindingSite> __Binding;

        DECLARE_TYPE( __Range,  Range);
        DECLARE_TYPE( __BindingSite, BindingSite );
        DECLARE_TYPE( __Binding, Binding);

    protected:
        // This function ensures that aComplexSpecies is sorted molwise.  It also 
        // ensures the two pieces of class data, the indexToMolMap and the 
        // molToIndexRangeMap are setup.
        void setupDataStructuresForCalculation(ComplexSpeciesRef aComplexSpecies);

        // This function returns a permutation such that, when applied to the complex
        // species, leaves it in a mol sorted state.  
        // Ie mol[i].getType() <= mol[j].getType() for all i<=k.
        Permutation calculateMolSortingPermutationForComplex(ComplexSpeciesCref aComplexSpecies);

        // This is the daddy function of them all.  This function calculates a Permutation such 
        // that when applied to the complex, will leave the complex in a state such that
        // the complex output state is minimal amongst all potential complex output states.
        Permutation calculateMinimizingPermutationForComplex(ComplexSpeciesCref aComplexSpecies);

        // ???
        void maximallyExtendPermutation(PermutationRef aRefPermutation, ComplexSpeciesCref aComplexSpecies);

        PartialTokenList calculatePartialTokenListForPermutation(ComplexSpeciesCref anAP, PermutationCref aPerm); 

        // Returns true if at least one permutation in the set of permutations is incomplete.
        // Returns false if every one is a complete and proper permutation/function.
        bool checkIfSetContainsIncompletePermutations(const std::set<Permutation>& setOfPPs) const;

        // Returns true if aComplexSpecies is sorted molwise.  Ie returns true
        // if for all i, j mol[i].getMolType() < mol[j].getMolType().
        bool checkPlexIsSortedByMol(ComplexSpeciesCref aComplexSpecies) const;


//         // TODO write brute force function and generally a checker to make sure all this works.
//         // This function works on molToIndex
//         void 
//         produceAllValidPermutations( std::set<Permutation>& setOfPerms);
        
//         template <typename T>
//         void 
//         extend( std::vector<T>& main, const std::vector<T>& sub)
//         {
//             for(typename std::vector<T>::const_iterator iter = sub.begin();
//                 iter != sub.end();
//                 ++iter)
//             {
//                 main.push_back( *iter);
//             }
//         }

//         template <typename T>
//         void
//         offset(std::set< std::vector<T> >& theSet, T aT)
//         {
//             for(typename std::set< std::vector<T> >::iterator iter = theSet.begin();
//                 iter != theSet.end();
//                 ++iter)
//             {
//                 offset( *iter, aT);
//             }
//         }

//         template <typename T>
//         void 
//         offset(typename std::vector<T>& main, T aT)
//         {
//             for(unsigned int i = 0;
//                 i != main.size();
//                 ++i)
//             {
//                 main[i] += aT;
//             }
//         }

//         void 
//         generate_Sn(std::set< std::vector<int> >& setOfPermutations, unsigned int i);
        

        std::map<int, std::string> indexToMolMap;
        std::map<std::string, Range > molToIndexRangeMap;

        struct MolIndexLessThanCmp : public std::binary_function<int, int, bool>
        {
            DECLARE_TYPE(ComplexSpecies::MolList, MolList);

            MolIndexLessThanCmp(ComplexSpeciesCref aComplexSpeciesForCmp);
            bool operator()(int ndx1, int ndx2);

        protected:
            MolListCref theComparisonMolList;
        };  

        struct namerBindingCmp : public std::binary_function<Binding, Binding, bool>
        {
            namerBindingCmp(int theBinding);
            bool operator()(Binding a, Binding b);

        protected:
            int bindingNumber;
        };

    };


    class BindingSorterByNdx
    {
    public:
        BindingSorterByNdx( unsigned int i):
            ndx(i)
        {}
        
        bool
        operator()(const ComplexSpecies::Binding& bindingTheFirst, const ComplexSpecies::Binding& bindingTheSecond) const
        {
            const ComplexSpecies::HalfBinding *firstHalfBinding, *secondHalfBinding;

            bindingTheFirst.first.first == ndx? firstHalfBinding = &bindingTheFirst.second : firstHalfBinding = &bindingTheFirst.first;
            

            bindingTheSecond.first.first == ndx? secondHalfBinding = &bindingTheSecond.second : secondHalfBinding = &bindingTheSecond.first;

            return *firstHalfBinding < *secondHalfBinding;
            
        }

    private:
        int ndx;
    };


} // namespace nmr



#endif 

