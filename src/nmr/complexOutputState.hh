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

#ifndef __COMPLEXOUTPUTSTATE_HH
#define __COMPLEXOUTPUTSTATE_HH

#include <vector>
#include <string>

namespace nmr
{

  namespace detail
  {

    struct ComplexOutputState
    {
      typedef std::string MolTokenStr;
      typedef std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > BindingTokenStr;
      typedef std::pair<std::string, std::pair<std::string, std::string> > ModificationTokenStr; 

      // ModificationTokens correspond 1-1 with modifications.  A ModToken 
      // (a, (b,c)) will have a and b as stringified ints.  In this case, a is the MolIndex,
      // b is the value of getModificationSiteIndex() and c in the value of the Modification.

      std::vector<MolTokenStr> theMolTokens;
      std::vector<BindingTokenStr> theBindingTokens;
      std::vector<ModificationTokenStr> theModificationTokens;
      
      void addMolTokenToOutputState(MolTokenStr aMolToken)
      {
	theMolTokens.push_back(aMolToken);
      }
      void addBindingTokenToOutputState(BindingTokenStr aBindingToken)
      {
	theBindingTokens.push_back(aBindingToken);
      }

      void addModificationTokenToOutputState(ModificationTokenStr aModificationToken)
      {
	theModificationTokens.push_back(aModificationToken);
      }

      void clear()
      {
	theMolTokens.clear();
	theBindingTokens.clear();
	theModificationTokens.clear();
      }
      
    };

  }

}



#endif
