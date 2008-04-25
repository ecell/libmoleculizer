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
    DECLARE_CLASS(ComplexSpeciesOutputMinimizer);
    class ComplexSpeciesOutputMinimizer
    {

    public:

        ComplexOutputState getMinimalOutputState(ComplexSpeciesCref aComplexSpecies);

    protected:
        typedef std::pair<int, int> __Range;
        DECLARE_TYPE( __Range,  Range);

        typedef std::pair<int, int> __BindingSite;
        DECLARE_TYPE( __BindingSite, BindingSite );

        typedef std::pair<BindingSite, BindingSite> __Binding;
        DECLARE_TYPE( __Binding, Binding);

        typedef std::set<Permutation>::const_iterator ConstPermSetIter;

        Permutation getPlexSortingPermutation(ComplexSpecies aComplexSpecies);
        Permutation calculateMinimizingPermutationForComplex(ComplexSpecies aComplexSpecies);
        void setupDataStructuresForCalculation();
        void maximallyExtendPermutation(Permutation& aRefPermutation);
        PartialTokenList calculatePartialTokenListForPermutation(ComplexSpecies& anAP, Permutation& aPerm); 
        Permutation calculateMolSortingPermutationForComplex(ComplexSpecies& aComplexSpecies);
        bool checkMolsInComplexAreIsSorted();
        bool checkExistsIncompletePermutations(std::set<Permutation>& setOfPPs);
        bool checkPlexIsSorted(const ComplexSpecies& aComplexSpecies) const;

        ComplexSpecies theUnnamedComplex;
        std::map<int, std::string> indexToMolMap;
        std::map<std::string, Range > molToIndexRangeMap;


        struct MolIndexLessThanCmp : public std::binary_function<int, int, bool>
        {
            DECLARE_TYPE(ComplexSpecies::MolList, MolList);

            MolIndexLessThanCmp(ComplexSpeciesCref aComplexSpeciesForCmp)
                : 
                theComparisonMolList( aComplexSpeciesForCmp.getMolList() )
            {}

            bool operator()(int ndx1, int ndx2)
            {
                return *theComparisonMolList[ndx1] < *theComparisonMolList[ndx2];
            } 
    

        protected:
            MolListCref theComparisonMolList;
        }; 

        struct namerBindingCmp : public std::binary_function<Binding, Binding, bool>
        {
            namerBindingCmp(int theBinding) 
                : 
                bindingNumber(theBinding)
            {}

            bool operator()(Binding a, Binding b)
            {
                // TODO/ 10 ???
                //This will sort the bindings by the binding number of a particular binding.  
	
                int aNumToComp;
                int bNumToComp;
	
                if (a.first.first==bindingNumber)
                {
                    aNumToComp=a.first.second;
                }
                else 
                {
                    aNumToComp=a.second.second;
                }
	
                if (b.first.first==bindingNumber)
                {
                    bNumToComp=b.first.second;
                }
                else 
                {
                    bNumToComp=b.second.second;
                }

                return (aNumToComp < bNumToComp);
            }
            int bindingNumber;

        };

    };

} // namespace nmr



#endif 

