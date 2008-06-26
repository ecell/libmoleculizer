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

#include <iterator>
#include <set>
#include <vector>
#include <utility>
#include <string>
#include <boost/foreach.hpp>
#include "complexSpeciesOutputMinimizer.hh"

namespace nmr
{

    ComplexOutputState 
    ComplexSpeciesOutputMinimizer::getMinimalOutputState(ComplexSpecies theComplexSpecies)
    { 

        // This ensures that theComplexSpecies is sorted and that the data structures, 
        // indexToMolMap and molToIndexRangeMap, are set up properly.
        setupDataStructuresForCalculation( theComplexSpecies );

        // This here is the money shot, where all the action occurs.  
        Permutation theUnnamedComplexMinimalOutputStatePermutation = 
            calculateMinimizingPermutationForComplex(theComplexSpecies);

        theComplexSpecies.applyPermutationToComplex(theUnnamedComplexMinimalOutputStatePermutation);

        ComplexOutputState theCanonicalOutputState;
        theComplexSpecies.constructOutputState(theCanonicalOutputState);

        return theCanonicalOutputState;
    }


    ComplexOutputState
    ComplexSpeciesOutputMinimizer::getMinimalOutputStateViaBruteForce( ComplexSpecies theComplexSpecies)
    {
        setupDataStructuresForCalculation( theComplexSpecies);
        std::set<Permutation> allPermutations;

        std::vector<unsigned int> signature;

        typedef const std::pair<std::string, Range> IteratedType;

        BOOST_FOREACH(IteratedType& i, molToIndexRangeMap )
        {
            signature.push_back( i.second.second - i.second.first);
        }

        bool DEBUG = false;
        if (signature[0] == 2 && signature[1] == 1 && signature[2] == 1)
        {
            DEBUG = true;
        }

        Permutation::generateAllPermutationsMatchingSignature( allPermutations, signature);



        PartialTokenList leastPartialTokenList = 
            calculatePartialTokenListForPermutation(theComplexSpecies, *allPermutations.begin());

        Permutation leastPermutation = *allPermutations.begin();



        for(std::set<Permutation>::const_iterator iter = allPermutations.begin();
            iter != allPermutations.end();
            ++iter)
        {
            PartialTokenList tmpPartialTokenList = 
                calculatePartialTokenListForPermutation( theComplexSpecies, *iter);
            if(tmpPartialTokenList < leastPartialTokenList)
            {
                leastPartialTokenList = tmpPartialTokenList;
                leastPermutation = *iter;
            }
        }

        theComplexSpecies.applyPermutationToComplex( leastPermutation );
        ComplexOutputState theCanonicalOutputState;

        theComplexSpecies.constructOutputState( theCanonicalOutputState );
        
        return theCanonicalOutputState;
    }

    Permutation 
    ComplexSpeciesOutputMinimizer::calculateMinimizingPermutationForComplex(ComplexSpeciesCref aComplexSpecies)
    {
        if( !checkPlexIsSortedByMol(aComplexSpecies) ) 
            throw nmr::GeneralNmrXcpt("Error in ComplexSpeciesOutputMinimizer::calculateMinimizingPermutationForComplex():\n The precondition that the ComplexSpecies be sorted is not met.");


        std::set<Permutation> setOfPossibleSolutions; // This is a set of partial permutations, which represent, all possible
                                                      // permutations (each incomplete pp represents a particular coset of
                                                      // permutations)

        setOfPossibleSolutions.insert( Permutation( aComplexSpecies.getNumberOfMolsInComplex()) );


        std::vector<PermutationName> theNextIterationsPartialPermutations; // Temorary storage for holding permutations which will
                                                                           // be in the setOfPossibleSolutions set during the next
                                                                           // iteration.


        //Basic idea of this algorithm as follows:
        //
        // 1.  For each permutation in setOfPossibleSolutions, figure out the least int not in that permutation.
        // 2.  That int will be a associated with the Mol that will ultimatly map to it.  Use this information to find all the
        //     indexes that can map to that least value.
        // 3.  Now prune the index list above (we are trying to find the smallest set of indexes that can extend the
        //     permutation in 1.) by only keeping indexes which correspond to a mol that has the lowest binding amonst all the
        //     corresponding mols to the cumulative indexes.
        // 4.  Now generate all the permutations that come from combining the permutation in question with each one of it's
        //     possible continuations.
        // 5.  For each one of these new permutations, extend it maximally.
        // 6.  Now generate the name corresponding to each permutation.  Find the least and throw out any permutation whose
        //     name isn't one of the least.
        // 7.  Repeat this process until there is no permutation in the set of permutations that has any negative entry.


        while( checkIfSetContainsIncompletePermutations(setOfPossibleSolutions) ) // 7.
        {
            theNextIterationsPartialPermutations.clear();

            //For each partial permutation, figure out how to automatically extend it into the next round.
            for(std::set<Permutation>::const_iterator partialPermIter = setOfPossibleSolutions.begin();
                partialPermIter != setOfPossibleSolutions.end();
                ++partialPermIter)
            {   
                //For each permutation in setOfPossibleSolutions, determine if it is complete or not.
                //If it is not, find all continuations, construct and extend each one, and add it to theNextIterationsPartialPermutations;
                //If it is, just add it as is to the partials of this round.  

                if ( !partialPermIter->getIsComplete() )
                {
                    const Permutation& refPerm = *partialPermIter;

                    //1. Find the most significant int that hasn't been fixed already.
                    int leastIntNotInPartial = partialPermIter->getLeastValueNotInPermutation();
		
                    //2. Find all the free mols of the associated type and put them in validIndexes
                    std::string molTypeToBeAdded = indexToMolMap[leastIntNotInPartial];
                    Range molRange = molToIndexRangeMap[molTypeToBeAdded];

                    std::vector<int> validIndexes;  
                    for(int i = molRange.first; i != molRange.second; ++i)
                    {
                        if ( partialPermIter->getValueAtPosition(i) == Permutation::UNDEF )
                        {
                            validIndexes.push_back(i);
                        }
                    }

                    //3.5
                    for(std::vector<int>::iterator i = validIndexes.begin();
                        i!=validIndexes.end();
                        ++i)
                    {
                        PermutationName tmpPn;
                        Permutation tmpPermutation(*partialPermIter, *i, leastIntNotInPartial);

                        //3.7
                        maximallyExtendPermutation(tmpPn.thePermutation, aComplexSpecies);
                        theNextIterationsPartialPermutations.push_back(tmpPn);

                    }
                }
                else
                {
                    // If the permutation is complete, just copy it...
                    PermutationName tmp;
                    tmp.thePermutation = *partialPermIter;
                    theNextIterationsPartialPermutations.push_back(tmp);
                }

            } //ends for partialPermIter loop, The creation of new partials for the round is complete


            //We begin the processing of the new partials by clearing setOfPossibleSolutions
            setOfPossibleSolutions.clear();
      
            //Generate the name for each Permutation
            for(std::vector<PermutationName>::iterator iter = theNextIterationsPartialPermutations.begin();
                iter != theNextIterationsPartialPermutations.end();
                ++iter)
            {
                iter->theCorrespondingPartialTokenList = calculatePartialTokenListForPermutation(aComplexSpecies, iter->thePermutation);
            }

            std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
            std::cout<< "Sorting " << theNextIterationsPartialPermutations.size() << std::endl;


             // int i = 0;
//              BOOST_FOREACH(const PermutationName& pm , theNextIterationsPartialPermutations)
//              {
//                  std::cout << i++<< ":" << std::endl;
//                  pm.print();
//                  std::cout << pm.theCorrespondingPartialTokenList << std::endl;
                
//              }
             std::cout << "###################################################################" << std::endl;
//           std::cout << "***" << std::endl;



            std::sort(theNextIterationsPartialPermutations.begin(),
                      theNextIterationsPartialPermutations.end());

            PartialTokenList leastPartial = theNextIterationsPartialPermutations.begin()->theCorrespondingPartialTokenList;
            for(std::vector<PermutationName >::iterator i = theNextIterationsPartialPermutations.begin();
                i != theNextIterationsPartialPermutations.end();
                ++i)
            {
                if (i->theCorrespondingPartialTokenList == leastPartial)
                {
                    setOfPossibleSolutions.insert( i->thePermutation );
                }
                else
                {
                    break;
                }
            }
        } //end while loop 

        return *(setOfPossibleSolutions.begin());
    }


    void ComplexSpeciesOutputMinimizer::setupDataStructuresForCalculation(ComplexSpeciesRef aComplexSpecies)
    {

        // First sort the UnnamedComplex by mol-types.  
        Permutation aMolSortingPermutation 
            = calculateMolSortingPermutationForComplex(aComplexSpecies);
        aComplexSpecies.applyPermutationToComplex(aMolSortingPermutation);

        // Clear out things.
        indexToMolMap.clear();
        molToIndexRangeMap.clear();

        //Initialize the indexToMolMap.
        std::vector<std::string> sortedMolNames;
        ComplexSpecies::MolList aMolList = aComplexSpecies.getMolList();

        for(unsigned int i = 0; 
            i!=aMolList.size(); 
            ++i)
        {
            sortedMolNames.push_back( aMolList[i]->getMolType() );
            indexToMolMap.insert(std::make_pair(i, aMolList[i]->getMolType()));
        }


        std::set<std::string> uniqueMolTypes(sortedMolNames.begin(),
                                             sortedMolNames.end());

        for(std::set<std::string>::iterator i = uniqueMolTypes.begin(); 
            i != uniqueMolTypes.end(); 
            ++i)
        {

            std::string currentMolType(*i);
            int first=0;
            int last;

            

            // TODO/5 Re-write this in a better fashion.  This is pretty ugly code.
            //Properly set the first index.
            while( !((first==0 && (sortedMolNames[first]==currentMolType)) || (first!=0 && sortedMolNames[first-1]!=currentMolType && sortedMolNames[first]==currentMolType)))
            {
                first++;
            }

            //Now set the last index.
            last=first+1;
            while(!(last==static_cast<int>(sortedMolNames.size()) || (sortedMolNames[last]!=sortedMolNames[last-1])))
            {
                last++;
            }
      
            molToIndexRangeMap.insert( std::make_pair( currentMolType, std::make_pair(first,last)));
        }

    }


    void ComplexSpeciesOutputMinimizer::maximallyExtendPermutation(Permutation& refPermToExtend, ComplexSpeciesCref aComplexSpecies)
    {
        uint iii = 0;

        // We iterate over the inversePerm in order to ensure that things are done from most significant to least significant in terms of 
        // the naming order.  

         Permutation inversePerm = refPermToExtend.invertPermutation();
         for(unsigned int indexToExtend = 0;
             indexToExtend != inversePerm.getDimension();
             ++indexToExtend)
         {
             
             inversePerm.getValueAtPosition( indexToExtend );
             

             // Nothing has been assigned to this index yet.  However, if this is the last "slot" a particular molType can end
             // up in, we can find that mol and assign it to this location.
//              if (inversePerm.getValueAtPosition( indexToExtend ) == Permutation::UNDEF )
//              {
                 
                  // Because the mols in the pre-arranged and post-arranged forms always have the same index (ie the mol at 
                  // index N will always be of the same type in aComplexSpecies and aComplexSpecies.applyPermutationToComplex(thePerm)
                  // we can find the set of appropriate indexes and see if this is the last one.  
//                  const std::string& theMolType( aComplexSpecies.getMolList()[ indexToExtend ]->getMolType() );
//                  if ( molToIndexRangeMap[ theMolType ].second  == indexToExtend + 1 )
//                  {
                      // AHA!!! It is the last one we need to set.  So now we figure out which in the range
//                      for(unsigned int originalNdx = molToIndexRangeMap[ theMolType ].first;
//                          originalNdx != molToIndexRangeMap[ theMolType ].second;
//                          ++originalNdx)
//                      {
//                          if ( refPermToExtend[originalNdx] == Permutation::UNDEF )
//                          {
//                              refPermToExtend.setValueAtPosition(originalNdx, indexToExtend);
//                              inversePerm.setValueAtPosition(indexToExtend, originalNdx);
//                              ++iii;
//                              break;
//                          }
//                  }
//              }
              

             
             // At first glance, it appears that this if block could be turned into a else block for the previous 
             // if block.  However, it may be that in the previous if block a value is set, in which case this block
             // should execute here as well.  Don't combine them.  
             if( inversePerm.getValueAtPosition(indexToExtend) != Permutation::UNDEF)
             {

                 int originalMolNdx( inversePerm.getValueAtPosition( indexToExtend ) );
                 std::vector<ComplexSpecies::Binding> theAppropriateBindings;


                 BOOST_FOREACH( const ComplexSpecies::Binding& theBinding, aComplexSpecies.getBindingList() )
                 {
                     if ( theBinding.first.first == indexToExtend || theBinding.second.first == indexToExtend ) 
                     { 
                         theAppropriateBindings.push_back( theBinding );
                     }
                 }

                 // Sort 'theAppropriateBindings' by the index of the non-participating binding...
                 std::sort( theAppropriateBindings.begin(),
                            theAppropriateBindings.end(),
                            BindingSorterByNdx( originalMolNdx ) );

                 BOOST_FOREACH( const ComplexSpecies::Binding& theBinding, theAppropriateBindings)
                 {
                     int molBindingPartnerNdx;
                     theBinding.first.first == originalMolNdx? 
                         molBindingPartnerNdx = theBinding.second.first :
                         molBindingPartnerNdx = theBinding.first.first;
                     
                     if (refPermToExtend.getValueAtPosition( molBindingPartnerNdx ) == Permutation::UNDEF )
                     {
                         // This is a mol whose value in the 'refPermToExtend' we can set.
                         
                         std::string theMolType( aComplexSpecies.getMolList()[molBindingPartnerNdx]->getMolType() );
                         int presumptiveIndex = molToIndexRangeMap[ theMolType ].first;
                         while(1)
                         {
                             if (inversePerm.getValueAtPosition( presumptiveIndex ) == Permutation::UNDEF)
                             {
                                 refPermToExtend.setValueAtPosition(molBindingPartnerNdx, presumptiveIndex);
                                 inversePerm.setValueAtPosition(presumptiveIndex, molBindingPartnerNdx);
                                 ++iii;
                                 break;
                                 
                             }
                             ++presumptiveIndex;
                         }
                         

                     }
                     
                 }
             }
         }
         

         std::cout << "In maximallyExtendPermutation " << iii << " indexes were extended..." << std::endl;
         
    }
    
    
             
    
    
         
             
                    

    PartialTokenList 
    ComplexSpeciesOutputMinimizer::calculatePartialTokenListForPermutation(ComplexSpeciesCref aComplexSpecies, PermutationCref aPerm)
    {

        PartialTokenList aComplexSpeciesPartialTokenList;
        aComplexSpeciesPartialTokenList.isComplete = true;

        // The mols are already sorted in an abstractPlex, so we can just copy those into aComplexSpeciesPartialTokenList.theMols.
        ComplexSpecies::MolList aComplexSpeciesMolList = aComplexSpecies.getMolList();

        copy(aComplexSpeciesMolList.begin(),
             aComplexSpeciesMolList.end(),
             back_inserter(aComplexSpeciesPartialTokenList.theMols));


        std::vector<Binding> aComplexSpeciesPartialTokenListBindings;

        typedef int MolNdx;
        typedef std::string BindingSite;
        typedef int BndNdx;
        typedef std::pair<std::pair<MolNdx, BndNdx>, std::pair<MolNdx, BndNdx> > Binding;
        typedef std::vector<Binding> BindingList;
      
        const BindingList& theBindings= aComplexSpecies.getBindingList();

        for(BindingList::const_iterator i = theBindings.begin();
            i != theBindings.end();
            ++i)
        {
            int f_index, f_bndNum, s_index, s_bndNum;
            std::string f_bndState, s_bndState;
      
            f_index = (*i).first.first;
            f_bndNum = (*i).first.second;
      
            s_index = (*i).second.first;
            s_bndNum = (*i).second.second;


            f_index = aPerm.getValueAtPosition(f_index);
            s_index = aPerm.getValueAtPosition(s_index);
            if (!(f_index == -1) && !(s_index == -1))
            {
                //then this entire binding is valid.
                Binding tmpBinding;
                if(f_index < s_index)
                {
                    tmpBinding.first.first = f_index;
                    tmpBinding.first.second = f_bndNum;


                    tmpBinding.second.first = s_index;
                    tmpBinding.second.second = s_bndNum;

                }
                else
                {
                    tmpBinding.first.first = s_index;
                    tmpBinding.first.second = s_bndNum;

                    tmpBinding.second.first = f_index;
                    tmpBinding.second.second = f_bndNum;

                }

                aComplexSpeciesPartialTokenListBindings.push_back(tmpBinding);
            }
            else 
            {
                aComplexSpeciesPartialTokenList.isComplete = false; 
            }

        }

        std::sort(aComplexSpeciesPartialTokenListBindings.begin(),
                  aComplexSpeciesPartialTokenListBindings.end());
  
        std::copy(aComplexSpeciesPartialTokenListBindings.begin(),
                  aComplexSpeciesPartialTokenListBindings.end(),
                  std::back_inserter(aComplexSpeciesPartialTokenList.theBindings));

        if (!aComplexSpeciesPartialTokenList.isComplete) 
        {
            return aComplexSpeciesPartialTokenList;
        }

        // I'm  not sure that this is stricly necessary, but at this moment, I'm worrying that it might be possible for an incomplete permutation to perfectly 
        // describe all the bindings and yet not describe all the modifications.  So, playing it safe, this may be more complicated than it needs to be.

    
        std::vector<std::pair<int, std::pair<std::string, std::string> > > tmpModList;
        Permutation inversePerm = aPerm.invertPermutation();



        const ComplexSpecies::MolList& theMols = aComplexSpecies.getMolList();
        for(unsigned int molNdx = 0; 
            molNdx != theMols.size(); 
            ++molNdx)
        {
            Mol::ModificationList currentModNdxMolList = theMols[molNdx]->getModificationList();

            for(Mol::ModificationList::iterator i = currentModNdxMolList.begin();
                i != currentModNdxMolList.end();
                ++i)
            {
                std::pair<int, std::pair<std::string, std::string> > unambiguousModification;
                unambiguousModification.first = inversePerm[molNdx];
                unambiguousModification.second = *i;
                tmpModList.push_back(unambiguousModification);
            }
	    
        }

        std::sort(tmpModList.begin(),
                  tmpModList.end());
        aComplexSpeciesPartialTokenList.theModifications.swap(tmpModList);

        return aComplexSpeciesPartialTokenList;
    }

    Permutation ComplexSpeciesOutputMinimizer::calculateMolSortingPermutationForComplex(ComplexSpeciesCref aComplexSpecies)
    {
        //initialize a vector with entries 0...size-1
        std::vector<int> permutationVect;
        for(unsigned int i=0; 
            i != aComplexSpecies.getNumberOfMolsInComplex(); 
            ++i)
        {
            permutationVect.push_back(i);
        }

        //Sort permutation using a function which compares the theMols belonging at that index.
        std::sort(permutationVect.begin(),
                  permutationVect.end(),
                  MolIndexLessThanCmp(aComplexSpecies) );

        Permutation inversePerm = Permutation(permutationVect);
        Permutation forwardPerm = inversePerm.invertPermutation();
        return forwardPerm;
    }

    bool ComplexSpeciesOutputMinimizer::checkIfSetContainsIncompletePermutations( const std::set<Permutation>& setOfPPs) const
    {
        // Iterate over the set and if any member is incomplete, 
        // return true; otherwise, return false.

        const bool SET_CONTAINS_INCOMPLETE_ELEMENT( true );
        const bool SET_DOES_NOT_CONTAIN_INCOMPLETE_ELEMENT( false );

        for(std::set<Permutation>::const_iterator iter = setOfPPs.begin(); 
            iter != setOfPPs.end(); 
            ++iter)
        {
            if (iter->getIsIncomplete()) return SET_CONTAINS_INCOMPLETE_ELEMENT;
        }

        return SET_DOES_NOT_CONTAIN_INCOMPLETE_ELEMENT;
    }


    bool 
    ComplexSpeciesOutputMinimizer::checkPlexIsSortedByMol(ComplexSpeciesCref aComplexSpecies) const
    {
        // This is the first time I've tried something like this.  It's probably too 
        // pedantic for most functions, but might be a nice trick to have around.
        const bool PLEX_IS_SORTED( true );
        const bool PLEX_IS_NOT_SORTED( false );

        if (aComplexSpecies.getNumberOfMolsInComplex() == 0)
        {
            throw GeneralNmrXcpt("Error in ComplexSpeciesOutputMinimizer::checkPlexIsSortedByMol.  Unnamed plex has not been initialized");
        }

        ComplexSpecies::MolListCref aMolList = aComplexSpecies.getMolList();

        if( aMolList.size() == 1) return PLEX_IS_SORTED;

        // This works because if for every pair of elements m_i and m_{i+1}
        // in a list, m_i <= m_{i+1}, the list is sorted.
        for(ComplexSpecies::MolList::const_iterator iter = aMolList.begin(), 
                iterPlusOne = iter + 1;
            iterPlusOne != aMolList.end();
            ++iter, ++iterPlusOne)
        {
            if (((*iter)->getMolType())>((*iterPlusOne)->getMolType()))
            {
                return PLEX_IS_NOT_SORTED;
            }
        }

        return PLEX_IS_SORTED;
    }


    ComplexSpeciesOutputMinimizer::MolIndexLessThanCmp::MolIndexLessThanCmp( ComplexSpeciesCref aComplexSpeciesForCmp)
        : 
        theComparisonMolList( aComplexSpeciesForCmp.getMolList() )
    {}

    bool
    ComplexSpeciesOutputMinimizer::MolIndexLessThanCmp::operator()(int ndx1, int ndx2)
    {
        return *theComparisonMolList[ndx1] < *theComparisonMolList[ndx2];
    }
    
    ComplexSpeciesOutputMinimizer::namerBindingCmp::namerBindingCmp( int theBinding)
        :
        bindingNumber( theBinding)
    {}


    bool ComplexSpeciesOutputMinimizer::namerBindingCmp::operator()( Binding a, Binding b)
    {
        // TODO/ 10 Rewrite/document/figure out WTF this does/etc.
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

    
    

//     void 
//     ComplexSpeciesOutputMinimizer::produceAllValidPermutations( std::set<Permutation>& setOfPermutations)
//     {
//         // This is a debugging function.  Given that the function "setupDataStructuresForCalculation"
//         // has been called, the member molToIndexRangeMap will be in effect.  
        
        
//         std::set<std::vector<int> > initialSet;
//         int second = molToIndexRangeMap.begin()->second.second;
//         int first = molToIndexRangeMap.begin()->second.first;
//         int n = second - first;
//         generate_Sn( initialSet, n);

//         // If this is it, just copy the set into set of permutations.
//         if (molToIndexRangeMap.size() == 1)
//         {
//             for(std::set<std::vector<int> >::iterator i = initialSet.begin();
//                 i != initialSet.end();
//                 ++i)
//             {
//                 setOfPermutations.insert( Permutation( *i) );
//             }
//             return;
//         }
        

//         std::map<std::string, Range>::const_iterator iter = molToIndexRangeMap.begin();
//         iter++;
//         for(;
//             iter != molToIndexRangeMap.end();
//             ++iter)
//         {
//             // Copy initialSet to tmpSet...
//             std::set<std::vector<int> > tmpSet(initialSet.begin(), initialSet.end());
//             initialSet.clear();

//             std::set<std::vector<int> > newSet;
//             generate_Sn( newSet, iter->second.second - iter->second.first);
            
//             for( std::set<std::vector<int> >::iterator ii = tmpSet.begin();
//                  ii != tmpSet.end();
//                  ++ii)
//             {
//                 for( std::set<std::vector<int> >::iterator jj = newSet.begin();
//                      jj != newSet.end();
//                      ++jj)
//                 {
//                     std::vector<int> thePerm( ii->begin(), ii->end() );
//                     std::vector<int> extendor( jj->begin(), jj->end() );
                    
//                     offset( extendor, (int) thePerm.size() );
//                     extend(thePerm, extendor);
//                     initialSet.insert( thePerm );

//                 }
//             }

//         }

//         for(std::set<std::vector<int> >::const_iterator ii = initialSet.begin();
//             ii != initialSet.end();
//             ++ii)
//         {
//             setOfPermutations.insert( Permutation( *ii));            
//         }
        
//     }


//     void 
//     ComplexSpeciesOutputMinimizer::generate_Sn( std::set< std::vector<int> >& setOfPermutations, unsigned int n)
//     {
        
//         if (n == 0) throw 666;
//         if (n == 1)
//         {
//             std::vector<int> permutation;
//             permutation.push_back( 0 );
//             setOfPermutations.insert( permutation );
//             return;
//         }

//         generate_Sn( setOfPermutations, n - 1);

//         std::set< std::list<int> > setOfPermLists;
//         std::set< std::list<int> > nextSetOfPermLists;
//         for( std::set<std::vector<int> >::iterator iter = setOfPermutations.begin();
//              iter != setOfPermutations.end();
//              ++iter)
//         {
//             setOfPermLists.insert( std::list<int>( iter->begin(), iter->end() ) );
//         }

//         setOfPermutations.clear();

//         for(std::set<std::list<int> >::iterator iter = setOfPermLists.begin();
//             iter != setOfPermLists.end();
//             ++iter)
//         {
//             int numToExtendWith = iter->size();
//             std::list<int> newList( iter->begin(), 
//                                     iter->end());

//             std::list<int>::iterator newIter = newList.begin();

//             do
//             {
//                 newList.insert(newIter, numToExtendWith);

//                 setOfPermutations.insert( std::vector<int>( newList.begin(),
//                                                             newList.end() ));

//                 newList.erase( newIter );
//             }
//             while(newIter++ != newList.end() );

//         }
        
//     }

}



