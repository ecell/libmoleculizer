/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2007  Nathan Addy
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


namespace complexspecies
{

  template <class molT>
  ComplexSpecies<molT>::ComplexSpecies() 
    : 
    theMolAliasToNdxMap(),
    theMols(),
    theBindings()
  {
    ; // do nothing
  }

  template <class molT>
  ComplexSpecies<molT>::ComplexSpecies(const ComplexSpecies<molT>& aComplexSpecies) 
    : 
    theMolAliasToNdxMap(aComplexSpecies.theMolAliasToNdxMap.begin(), aComplexSpecies.theMolAliasToNdxMap.end()),
    theMols(aComplexSpecies.theMols.begin(), aComplexSpecies.theMols.end()),
    theBindings(aComplexSpecies.theBindings.begin(), aComplexSpecies.theBindings.end())
  {
    ; // do nothing
  }

  template <class molT>
  void ComplexSpecies<molT>::addMolToComplex(molT someMol, Alias anAlias)
  {
    MolMapIter loc= theMolAliasToNdxMap.find(anAlias);
    if (loc != theMolAliasToNdxMap.end()) 
      {
	std::string anErrorMessage = "Error: complexSpecies<molT>::addMolToComplex failed.\n" + anAlias + " already used in this complex.";
	throw CSXcpt(anErrorMessage);
      }

    theMols.push_back(someMol);

    int newMolIndex = theMols.size()-1;
    theMolAliasToNdxMap.insert( std::make_pair(anAlias, newMolIndex) );
  }

  template <class molT>
  void ComplexSpecies<molT>::addBindingToComplex(Alias firstMolAlias, 
						 BindingSite firstMolBindingSiteAlias, 
						 Alias secondMolAlias, 
						 BindingSite secondMolBindingSiteAlias)
  {
    MolMapIter locFirstMol = theMolAliasToNdxMap.find(firstMolAlias);
    MolMapIter locSecondMol = theMolAliasToNdxMap.find(secondMolAlias);

    if (locFirstMol == theMolAliasToNdxMap.end()) 
      {
	std::string anErrorMessage = "Error: complexSpecies<molT>::addBindingToComplex failed.\n" + firstMolAlias + " does not exist in this complex.";
	throw CSXcpt(anErrorMessage);
      }
    
    if (locSecondMol == theMolAliasToNdxMap.end()) 
      {
	std::string anErrorMessage = "Error: complexSpecies<molT>::addBindingToComplex failed.\n" + secondMolAlias + " does not exist in this complex.";
	throw CSXcpt(anErrorMessage);
      }

    MolNdx firstMolNdx = theMolAliasToNdxMap[firstMolAlias];
    MolNdx secondMolNdx = theMolAliasToNdxMap[secondMolAlias];

    if( !(theMols[firstMolNdx].checkIfBindingSiteExists(firstMolBindingSiteAlias)))
      {
	// throw exception
	std::string anErrorMessage = "Error: complexSpecies<molT>::addBindingToComplex failed.\nMol with alias " + firstMolAlias + " does not have binding site " + firstMolBindingSiteAlias;
	throw CSXcpt(anErrorMessage);
      }

    if( !(theMols[secondMolNdx].checkIfBindingSiteExists(secondMolBindingSiteAlias)))
      {
	// throw exception
	std::string anErrorMessage = "Error: complexSpecies<molT>::addBindingToComplex failed.\nMol with alias " + secondMolAlias + "does not have binding site " + secondMolBindingSiteAlias;
	throw CSXcpt(anErrorMessage);
      }

    if( theMols[firstMolNdx].checkIfBindingSiteIsBound(firstMolBindingSiteAlias))
      {
	// throw exception
	std::string anErrorMessage = "Error: complexSpecies<molT>::addBindingToComplex failed.\n Mol " + firstMolAlias + "is already bound at site " + firstMolBindingSiteAlias;
	throw CSXcpt(anErrorMessage);
      }

    if( theMols[secondMolNdx].checkIfBindingSiteIsBound(secondMolBindingSiteAlias))
      {
	// throw exception
	std::string anErrorMessage = "Error: complexSpecies<molT>::addBindingToComplex failed.\n Mol " + secondMolAlias + "is already bound at site " + secondMolBindingSiteAlias;
	throw CSXcpt(anErrorMessage);
      }

    BndNdx firstMolBindingNdx = theMols[firstMolNdx].getBindingSiteInteger(firstMolBindingSiteAlias);
    BndNdx secondMolBindingNdx = theMols[secondMolNdx].getBindingSiteInteger(secondMolBindingSiteAlias);

    Binding aBinding;
    aBinding.first.first = firstMolNdx;
    aBinding.first.second = firstMolBindingNdx;
    aBinding.second.first = secondMolNdx;
    aBinding.second.second = secondMolBindingNdx;
    
    theBindings.push_back( aBinding);
  }


  template <class molT>
  int ComplexSpecies<molT>::getNumberOfMolsInComplex() const
  {
    return theMols.size();
  }

  template <class molT>
  int ComplexSpecies<molT>::getNumberOfBindingsInComplex() const
  {
    return theBindings.size();
  }


  template <class molT>
  const typename ComplexSpecies<molT>::MolList& ComplexSpecies<molT>::getMolList() const
  {
    return theMols;
  }

  template <class molT>
  const typename ComplexSpecies<molT>::BindingList& ComplexSpecies<molT>::getBindingList() const
  {
    return theBindings;
  }


  template <class molT>
  void ComplexSpecies<molT>::constructOutputState(detail::ComplexOutputState& anOutputState) const
  {
    anOutputState.clear();
  
    for(typename std::vector<molT>::const_iterator ii = theMols.begin();
	ii != theMols.end();
	++ii)
      {

	MolTokenStr aMolToken = ii->getMolType();
	anOutputState.addMolTokenToOutputState(aMolToken);
      }


    //Copy in the bindings.
    //The bindings are pairs of pair<int, int>, so we must unroll it ourselves.
    for(BindingList::const_iterator i=theBindings.begin();
	i!=theBindings.end();
	++i)
      {

	BindingTokenStr aBindingToken;
	
	int firstfirst, firstsecond, secondfirst, secondsecond;
	firstfirst = (*i).first.first;
	firstsecond = (*i).first.second;
	secondfirst = (*i).second.first;
	secondsecond = (*i).second.second;

	aBindingToken.first.first = detail::stringify(firstfirst);
	aBindingToken.first.second = detail::stringify(firstsecond);
	aBindingToken.second.first = detail::stringify(secondfirst);
	aBindingToken.second.second = detail::stringify(secondsecond);

	anOutputState.addBindingTokenToOutputState(aBindingToken);
      }

    for(int molNdx=0; molNdx != this->getNumberOfMolsInComplex(); ++molNdx)
      {

	std::string strMolNdx = detail::stringify(molNdx);
	ModificationList currentMolModificationList = theMols[molNdx].getModificationList();
	
	for(unsigned int modificationNdx = 0; 
	    modificationNdx != currentMolModificationList.size(); 
	    ++modificationNdx)
	  {
	    std::string strModificationSite = currentMolModificationList[modificationNdx].first;
	    std::string aModificationValue = currentMolModificationList[modificationNdx].second;
	    
	    ModificationTokenStr aModificationToken;

            // This should be stringified.

	    aModificationToken.first = detail::stringify(molNdx);


	    aModificationToken.second.first = strModificationSite;
	    aModificationToken.second.second = aModificationValue;

	    anOutputState.addModificationTokenToOutputState(aModificationToken);
	  }
      }
    
    return;
  }

  template <class molT>
  void ComplexSpecies<molT>::constructPartialTokenList(detail::PartialTokenList<molT>& rComplexPartialTokenList) const
  {
    for(typename MolList::const_iterator index = theMols.begin();
	index != theMols.end();
	++index)
      {
	rComplexPartialTokenList.theMols.push_back(*index);
      }
    for(typename MolList::const_iterator index = theBindings.begin();
	index != theBindings.end();
	++index)
      {
	rComplexPartialTokenList.theBindings.push_back(*index);
      }

    for(int molNdx = 0; molNdx != theMols.size(); ++molNdx)
      {
	ModificationList currentModNdxMolList = theMols[molNdx].getModificationList();

	for(typename ModificationList::iterator i = currentModNdxMolList.begin();
	    i != currentModNdxMolList.end();
	    ++i)
	  {
	    std::pair<int, std::pair<std::string, std::string> > unambiguousModification;
	    unambiguousModification.first = molNdx;
	    unambiguousModification.second = *i;
	    rComplexPartialTokenList.push_back(unambiguousModification);
	  }
	    
      }
    
  }


  template <class molT>
  void ComplexSpecies<molT>::applyPermutationToComplex(detail::Permutation& aPermutation)
  {
    if (aPermutation.getPermutationSize()==0)
      {
	throw CSXcpt("ComplexSpecies is empty");
      }
    if (!aPermutation.getIsComplete())
      {
	throw CSXcpt("Permutation is not complete");
      }
    if( aPermutation.getPermutationSize() != this->getNumberOfMolsInComplex() )
      {
	throw CSXcpt("Permutation in not the same size as our mol");
      }

    detail::Permutation theInversePerm=aPermutation.invertPermutation();

    MolList updatedMolList;
    BindingList updatedBindingList;

    for(unsigned int molNdx=0; molNdx!=theMols.size(); ++molNdx)
      {
	updatedMolList.push_back( theMols[theInversePerm[molNdx]] );
      }

    for(BindingList::iterator bndIter = theBindings.begin();
	bndIter != theBindings.end();
	++bndIter)
      {
	Binding copyOfBinding = (*bndIter);
	Binding updatedBinding;

	updatedBinding.first.first = aPermutation[copyOfBinding.first.first];
	updatedBinding.first.second = copyOfBinding.first.second;
	updatedBinding.second.first = aPermutation[copyOfBinding.second.first];
	updatedBinding.second.second = copyOfBinding.second.second;
	
	sortBinding(updatedBinding);      

	updatedBindingList.push_back(updatedBinding);
      }

    sort(updatedBindingList.begin(),
	 updatedBindingList.end());

    theMols.swap(updatedMolList);
    theBindings.swap(updatedBindingList);

    //Finally update the theMolAliasToNdxMap
    for(std::map<Alias, int>::iterator i=theMolAliasToNdxMap.begin();
	i!=theMolAliasToNdxMap.end();
	++i)
      {
	int originalIndex = (*i).second;
	(*i).second = aPermutation[originalIndex];
      }
  }


  template <class molT>
  void ComplexSpecies<molT>::sortBinding(Binding& aBinding)
  {
    if (aBinding.first > aBinding.second)
      {
	swap(aBinding.first, aBinding.second);
      }
  }



}
