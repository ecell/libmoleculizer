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

#ifndef __NAMEASSEMBLER_HH
#define __NAMEASSEMBLER_HH

#include "nmrExceptions.hh"

#include "complexSpecies.hh"
#include "complexOutputState.hh"
#include "complexSpeciesOutputMinimizer.hh"

#include <iterator>
#include <vector>
#include <string>
#include <utility>

namespace nmr
{

  template <typename molT>
  class NameAssembler
  {
  public:
    NameAssembler(const std::string& assemblerName)
      :
      name(assemblerName)
    {
    }

    virtual ~NameAssembler()
    {}

    std::string createCanonicalName(ComplexSpecies<molT> aComplexSpecies) const
    {
      detail::ComplexSpeciesOutputMinimizer<molT> canonicalNameGenerator;
      detail::ComplexOutputState nameSkeleton = canonicalNameGenerator.getMinimalOutputState( aComplexSpecies );
      return createNameFromOutputState( nameSkeleton );
    }

    std::string createName(const ComplexSpecies<molT>& aComplexSpecies) const
    {
      detail::ComplexOutputState anOutputState;
      aComplexSpecies.constructOutputState( anOutputState );
      return createNameFromOutputState( anOutputState );
    }

    const std::string& getName()
    {
      return name;
    }

  protected:

    const std::string name;

    // This function is intended as a debugish sort of a function, which throws an exception if a complex 
    // output state does not get encoded and then decoded in the same way

    static void assertEncodeDecodeAccuracy(NameAssembler<molT>* ptrNameAssembler) throw(encodeDecodeInconsistencyXcpt) 
    {
      detail::ComplexOutputState aComplexOutputState;

      detail::ComplexOutputState decodedEncodingOS = ptrNameAssembler->createOutputStateFromName( ptrNameAssembler->createNameFromOutputState( aComplexOutputState ) );
      if (! (decodedEncodingOS == aComplexOutputState)  )
        {
          throw encodeDecodeInconsistencyXcpt(ptrNameAssembler->getName());
        }
    }



    virtual std::string createNameFromOutputState( const detail::ComplexOutputState& aComplexOutputState) const = 0;
    virtual detail::ComplexOutputState createOutputStateFromName( const std::string& aComplexOutputState) const = 0;
    

  };

}

#endif
