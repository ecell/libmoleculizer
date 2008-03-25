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

#include <map>
#include <functional>

#include "complexSpecies.hh"
#include "nameAssembler.hh"
#include "permutation.hh"
#include "partialTokenList.hh"
#include "permutationName.hh"

namespace nmr
{

    namespace detail
    {
        template <class molT>
        class ComplexSpeciesOutputMinimizer
        {

        public:
            typedef molT Mol;
            typedef std::pair<int, int> Range;
            typedef std::pair<int, int> BindingSite;
            typedef std::pair<BindingSite, BindingSite> Binding;
            typedef std::set<Permutation>::const_iterator ConstPermSetIter;  

        public:

            ComplexSpeciesOutputMinimizer();
            ComplexOutputState getMinimalOutputState(ComplexSpecies<molT> aComplexSpecies);

        protected:

            Permutation getPlexSortingPermutation(ComplexSpecies<molT> aComplexSpecies);
            Permutation calculateMinimizingPermutationForComplex(ComplexSpecies<molT> aComplexSpecies);
            void setupDataStructuresForCalculation();
            void maximallyExtendPermutation(Permutation& aRefPermutation);
            PartialTokenList<molT> calculatePartialTokenListForPermutation(ComplexSpecies<molT>& anAP, Permutation& aPerm); 
            Permutation calculateMolSortingPermutationForComplex(ComplexSpecies<molT>& aComplexSpecies);
            bool checkMolsInComplexAreIsSorted();
            bool checkExistsIncompletePermutations(std::set<Permutation>& setOfPPs);
            bool checkPlexIsSorted(const ComplexSpecies<molT>& aComplexSpecies) const;

        protected:

            ComplexSpecies<molT> theUnnamedComplex;
            std::map<int, std::string> indexToMolMap;
            std::map<std::string, Range > molToIndexRangeMap;


            struct MolIndexLessThanCmp : public std::binary_function<int, int, bool>
            {
                typedef std::vector<molT> MolList;

                MolIndexLessThanCmp(const ComplexSpecies<molT>& aComplexSpeciesForCmp) 
                    : 
                    theComparisonMolList(aComplexSpeciesForCmp.getMolList())
                {
                    ; // do nothing
                }

                bool operator()(int ndx1, int ndx2)
                {
                    return theComparisonMolList[ndx1] < theComparisonMolList[ndx2];
                } 
	
            protected:
                const MolList& theComparisonMolList;
            };
        };


        struct namerBindingCmp
        {

            typedef std::pair<int, int> BindingSite;
            typedef std::pair<BindingSite, BindingSite> Binding;
     
            namerBindingCmp(int theBinding) 
                : 
                bindingNumber(theBinding)
            {
                ; // do nothing
            }

            bool operator()(Binding a, Binding b)
            {
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



    }

} //this r-brace ends namespace nmr; make sure it is on this side of the next #include line.  Otherwise trouble.  

#include "complexSpeciesOutputMinimizerImpl.hh"

#endif 

