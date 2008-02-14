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

#ifndef __BASICNAMEASSEMBLER_HH
#define __BASICNAMEASSEMBLER_HH

#include <string>
#include <vector>
#include "nameAssembler.hh"
#include "complexSpecies.hh"
#include "complexOutputState.hh"

namespace complexspecies
{

  template <typename molT>
  class basicNameAssembler : public NameAssembler<molT>
  {
  public:

    std::string createNameFromOutputState( const detail::ComplexOutputState& aCOS) const
    {
      std::string name("");
      for(std::vector<std::string>::const_iterator i = aCOS.theMolTokens.begin();
          i != aCOS.theMolTokens.end();
          ++i)
        {
          name += *i + ", ";
        }

      name = name.substr(0 , name.length() - 2);
      name += " -- ";
      
      for( std::vector< detail::ComplexOutputState::BindingTokenStr >::const_iterator i = aCOS.theBindingTokens.begin();
           i != aCOS.theBindingTokens.end();
           ++i)
        {
          name +=  (*i).first.first + (*i).first.second + (*i).second.first + (*i).second.second + ", ";
        }

      if (aCOS.theBindingTokens.size())
        {
          name = name.substr(0 , name.length() - 2);
        }

      name += " -- ";

      
      for( std::vector< detail::ComplexOutputState::ModificationTokenStr >::const_iterator iter = aCOS.theModificationTokens.begin();
           iter != aCOS.theModificationTokens.end();
           ++iter)
        {
          name += "( " + (*iter).first + ", " + (*iter).second.first + ", " + (*iter).second.second + "), ";
        }
      
      return name;
    }

  };

}
#endif
