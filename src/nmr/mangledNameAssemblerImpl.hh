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

#include "utl/utility.hh"
#include "nmr/nmrExceptions.hh"
#include "complexSpeciesOutputMinimizer.hh"
#include <iostream>
#include <utility>

namespace nmr
{
  template <typename molT>
  detail::ComplexOutputState MangledNameAssembler<molT>::createOutputStateFromName( const std::string& name) const
  {
    std::vector<std::string> tokens;
    utl::tokenize(name, tokens, "___");

    std::string theMolString( tokens[0] );
    std::string theBindingString(tokens[1] );
    std::string theModificationString( tokens[2] );

    std::vector<detail::ComplexOutputState::MolTokenStr> molTokens;
    std::vector<detail::ComplexOutputState::BindingTokenStr> bindingTokens;
    std::vector<detail::ComplexOutputState::ModificationTokenStr> modificationTokens;
      
    parseMolString(theMolString, molTokens);
    parseBindingString( theBindingString, bindingTokens);
    parseModificationString( theModificationString, modificationTokens);

    detail::ComplexOutputState newOutputState;
    newOutputState.theMolTokens = molTokens;
    newOutputState.theBindingTokens = bindingTokens;
    newOutputState.theModificationTokens = modificationTokens;

    std::cout << "#############################################" << std::endl;
    std::cout << name << std::endl << "->" << std::endl << newOutputState << std::endl;
    return newOutputState;
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

  template <typename molT>
  void MangledNameAssembler<molT>::parseMolString(const std::string& molString, std::vector<detail::ComplexOutputState::MolTokenStr>& molTokenVector) const
  {
    std::string::size_type index = 0;
    while(true)
      {
        if( index == molString.size() ) return;
        
        std::string::size_type lengthOfMolNameToken(0);
        
        // Index here is going to point to the LengthToken
        if ( molString[index] == '_' )
          {
            index+=1;
            std::string::size_type numberEnd = molString.find( "_", index);
            utl::from_string<std::string::size_type>(lengthOfMolNameToken, std::string(molString, index, numberEnd - index));
            index = numberEnd + 1;
          }
        else
          {
            utl::from_string<std::string::size_type>(lengthOfMolNameToken, std::string(molString, index, 1));
            index += 1;
          }

        // Post conditions:
        // 1.  Index should point to the first character of the name
        // 2.  length should be the appropriate length of the name.

        molTokenVector.push_back( std::string(molString, index, lengthOfMolNameToken) );
        index += lengthOfMolNameToken;

        return;
      }
    
    
  }

  template <typename molT>
  void MangledNameAssembler<molT>::parseBindingString(const std::string& bindingString, std::vector<detail::ComplexOutputState::BindingTokenStr>& bindingTokenVector ) const
  {
    
    std::vector<std::string> tmpVector;
    std::string::size_type currentBindingLength;
    
    std::string::size_type index = 0;
    
    while(true)
      {
        if ( index == bindingString.size() ) break;
        
        if ( bindingString[index] == '_' )
          {
            index += 1;
            //
            std::string::size_type endOfLengthToken = bindingString.find("_", index);
            utl::from_string<std::string::size_type>(currentBindingLength, std::string(bindingString, index, endOfLengthToken - index));
            index = endOfLengthToken + 1;
          }
        else
          {
            currentBindingLength = 1;
          }

        // Post conditions:
        // 1. Index should point to the first character of the name
        // 2. currentBindingLength should be the length of the next index number.

        tmpVector.push_back( std::string(bindingString, index, currentBindingLength) );
        index += currentBindingLength;
      }

    // Checking for an badly formed binding string
    // Error condition.
    if (tmpVector.size() % 4 != 0)
      {
        throw badBindingNameXcpt( bindingString);
      }

    // Package the sequence of integers into a sequence of bindingTokens
    for(unsigned int ii = 0;
        ii != tmpVector.size() / 4;
        ++ii)
      {
        bindingTokenVector.push_back( std::make_pair( std::make_pair(tmpVector[4*ii], tmpVector[4*ii + 1]), 
                                                      std::make_pair(tmpVector[4*ii+2], tmpVector[4*ii+3]) ) );
      }

    return;
            
  }

  template <typename molT>
  void MangledNameAssembler<molT>::parseModificationString(const std::string& modString, std::vector<detail::ComplexOutputState::ModificationTokenStr>& modificationTokenVector) const
  {
    // A mod string should be ( string, (string, string) )

    std::vector<std::string> tmpVector;
    std::string::size_type index(0);
    std::string::size_type numberOfCharactersInNextToken;
    
    while(true)
      {
        if (index == modString.size() ) break;

        // According to processModificationToken, each token begins with a _
        if (modString[index] != '_') throw badModificationNameXcpt( modString );

        ++index;

        if (modString[index] == '_')
          {
            // We have a multi-charecter length to read.

            index += 1;
            
            std::string::size_type endOfLengthString = modString.find("_", index);
            utl::from_string<std::string::size_type>(numberOfCharactersInNextToken, 
                                                     std::string(modString, 
                                                                 index, 
                                                                 endOfLengthString - index));
            index = endOfLengthString + 1;
          }
        else
          {
            utl::from_string<std::string::size_type>(numberOfCharactersInNextToken,
                                                     std::string(modString, 
                                                                 index,
                                                                 1));
            index += 1;
          }

        tmpVector.push_back( std::string(modString,
                                         index, 
                                         numberOfCharactersInNextToken) );

        index += numberOfCharactersInNextToken;

      }

    if( tmpVector.size() % 3 != 0)
      {
        std::cerr << "##############################" << std::endl;
        std::cerr << "Error in parsing modificaton string '" << modString << "'." << std::endl;
        std::cerr << "Which was parsed as: " << std::endl;
        for (unsigned int ii = 0;
             ii != tmpVector.size();
             ++ii)
          {
            std::cerr << ii << ":\t\"" << tmpVector[ii] << "\"" << std::endl;
          }
        throw badModificationNameXcpt( modString );
      }

    for( unsigned int ii = 0;
         ii != tmpVector.size() / 3;
         ++ii)
      {
        modificationTokenVector.push_back( std::make_pair( tmpVector[3*ii],
                                                           std::make_pair(tmpVector[3*ii + 1], tmpVector[3*ii + 2])));
      }


  }


}
