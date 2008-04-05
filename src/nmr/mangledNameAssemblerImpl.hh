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


#include "complexSpeciesOutputMinimizer.hh"
#include <iostream>

namespace nmr
{
  template <typename molT>
  detail::ComplexOutputState MangledNameAssembler<molT>::createOutputStateFromName( const std::string& name) const
  {
    std::string deliminator("___");

    std::string::size_type firstCharOfMols = name.find_first_not_of(deliminator, 0);
    std::string::size_type lastCharOfMols = name.find_first_of(deliminator, firstCharOfMols);
    std::string::size_type firstCharOfBindings = name.find_first_not_of(deliminator, lastCharOfMols);
    std::string::size_type lastCharOfBindings = name.find_first_of(deliminator, firstCharOfBindings);
    std::string::size_type firstCharOfModifications = name.find_first_not_of(deliminator, lastCharOfBindings);
    std::string::size_type lastCharOfModifications = name.size();


    std::string molString(name,
                          firstCharOfMols, 
                          lastCharOfMols);
    
    std::string bindingsString(name,
                               firstCharOfBindings,
                               lastCharOfBindings);
    
    std::string modificationsString(name,
                                    firstCharOfModifications);
    
    std::cout << molString << std::endl;
    std::cout << bindingsString << std::endl;
    std::cout << modificationsString << std::endl;
    
  }
  
  template <typename molT>
  std::string MangledNameAssembler<molT>::createNameFromOutputState( const detail::ComplexOutputState& aComplexSpeciesOutputState) const
  {
    std::string aComplexSpeciesName("");

    aComplexSpeciesName += "___";
    aComplexSpeciesName += constructMangledMolList(aComplexSpeciesOutputState);
    aComplexSpeciesName += "___";
    aComplexSpeciesName += constructMangledBindingList(aComplexSpeciesOutputState);
    aComplexSpeciesName += "___";
    aComplexSpeciesName += constructMangledModificationList(aComplexSpeciesOutputState);

    createOutputStateFromName( aComplexSpeciesName );

    return aComplexSpeciesName;
  }

  template <typename molT>
  std::string MangledNameAssembler<molT>::constructMangledMolList(const detail::ComplexOutputState& aComplexOutputState) const 
  {
    std::string mangledMolName("");
   
    std::vector<std::string>::const_iterator molIterator;
    for(molIterator = aComplexOutputState.theMolTokens.begin();
	molIterator != aComplexOutputState.theMolTokens.end();
	++molIterator)
      {
 	std::string tmpName(*molIterator);
	tmpName=getEncodedLength(tmpName)+tmpName;
	mangledMolName+=tmpName;
      }

    return mangledMolName;
  }


  template <typename molT>
  std::string MangledNameAssembler<molT>::constructMangledBindingList(const detail::ComplexOutputState& aComplexOutputState) const
  {
    std::string mangledBindingList("");

    std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > >::const_iterator bndIterator;
    for( bndIterator = aComplexOutputState.theBindingTokens.begin();
	 bndIterator != aComplexOutputState.theBindingTokens.end();
	++bndIterator)
      {
	std::string bindingStr;
	std::string first(bndIterator->first.first);
	std::string second(bndIterator->first.second);
	std::string third(bndIterator->second.first);
	std::string fourth(bndIterator->second.second);
	
	bindingStr = processBindingString(first) + processBindingString(second) + processBindingString(third) + processBindingString(fourth);
	mangledBindingList += bindingStr;
       }

     return mangledBindingList;
  }
  
  template <typename molT>
  std::string MangledNameAssembler<molT>::constructMangledModificationList(const detail::ComplexOutputState& aComplexOutputState) const
  {
     std::string mangledModString("");
     
     typedef std::pair<std::string, std::pair<std::string, std::string> > Modification;
     typedef std::vector<Modification> ModificationList;

     for(ModificationList::const_iterator i = aComplexOutputState.theModificationTokens.begin();
	 i!=aComplexOutputState.theModificationTokens.end();
	 ++i)
       {

	 std::string currentModificationString;
	 std::string first = i->first;
	 std::string second = i->second.first;
	 std::string third = i->second.second;
	 
	 currentModificationString = processModificationToken(first) + processModificationToken(second) + processModificationToken(third);
	 mangledModString += currentModificationString;
       }

    return mangledModString;
  }


  template <typename molT>
  std::string MangledNameAssembler<molT>::getEncodedLength(const std::string& stringInQuestion) const
  {
    int len=stringInQuestion.size();
    if (len<10)
      {
	//if len is between 1 and 9 inclusive then it is only one charector long
	return detail::stringify(len);
      }
    else
      {
	std::string tmp;
	tmp="_" + detail::stringify(len) + "_";
	return tmp;
      }
  }


}
