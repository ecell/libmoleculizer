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

#ifndef __READABLENAMEASSEMBLER_HH
#define __READABLENAMEASSEMBLER_HH

#include <string>
#include <vector>
#include "nameAssembler.hh"
#include "complexSpecies.hh"
#include "complexOutputState.hh"

namespace nmr
{
  template <typename molT>
  class readableNameAssembler : public NameAssembler<molT>
  {
    
  public:
    std::string createNameFromOutputState( const detail::ComplexOutputState& aCOS) const
    {
      // We would like the output to be something like the following.
      // X(phosphorylated):Y:Z(not_phos).

      std::string name("");

      for(unsigned int molNdx = 0; molNdx != aCOS.theMolTokens.size(); ++molNdx)
        {
          name += aCOS.theMolTokens[molNdx];

          std::string modSitesString("");
          int numberModSites = 0;
          std::string molNdxAsString = this->stringify( molNdx );
          
          for( unsigned int i = 0 ; i!= aCOS.theModificationTokens.size(); ++i)
            {
              if (aCOS.theModificationTokens[i].first == molNdxAsString)
                {
                  modSitesString += aCOS.theModificationTokens[i].second.second;
                  modSitesString += ",";
                  numberModSites++;
                }
            }
          
          if (numberModSites)
            {
              modSitesString = modSitesString.substr(0, modSitesString.length() - 1);
              modSitesString = "(" + modSitesString + ")";
              name += modSitesString;
            }
          name += "-";
        }
      name = name.substr(0, name.length() - 1);
      return name;
    }

  private:
    std::string 
    stringify(int i) const
    {
      std::stringstream stream;
      stream << i;
      
      std::string numAsString;
      stream >> numAsString;
      
      return numAsString;
    }

  };

}

#endif
