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

#include <iterator>
#include <vector>
#include <string>
#include <utility>

#include "complexSpecies.hh"
#include "complexOutputState.hh"
#include "complexSpeciesOutputMinimizer.hh"

namespace nmr
{

  template <typename molT>
  class NameAssembler
  {
  public:
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

    virtual std::string createNameFromOutputState( const detail::ComplexOutputState& aComplexOutputstate) const = 0;

  };

}

#endif
