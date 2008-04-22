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


#include "utl/macros.hh"

#include <vector>
#include <string>
#include <iostream>

namespace nmr
{

    // A Complex Output State is an intermediate between a complex species and a mangled
    // name 

    DECLARE_CLASS( ComplexOutputState );
    struct ComplexOutputState
    {
        DECLARE_TYPE( std::string, 
                      MolTokenStr);

        typedef std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > __BindingTokenStr;
        typedef std::pair<std::string, std::pair<std::string, std::string> > __ModificationTokenStr;

        DECLARE_TYPE( __BindingTokenStr, BindingTokenStr);
        DECLARE_TYPE( __ModificationTokenStr, ModificationTokenStr);

        std::vector<MolTokenStr> theMolTokens;
        std::vector<BindingTokenStr> theBindingTokens;

        // ModificationTokens correspond 1-1 with modifications.  A ModToken 
        // (a, (b,c)) will have a and b as stringified ints.  In this case, a is the MolIndex,
        // b is the value of getModificationSiteIndex() and c in the value of the Modification.
        std::vector<ModificationTokenStr> theModificationTokens;

        bool operator!=(const ComplexOutputState& other) const
        {
            return !( *this == other );
        }

        bool operator==(const ComplexOutputState& other) const
        {

//       // If any of the sizes of the three vectors of tokens don't match up, 
//       // they cannot be equal.
//       if (theMolTokens.size() != other.theMolTokens.size() || 
//           theBindingTokens.size() != other.theBindingTokens.size() || 
//           theModificationTokens.size() != other.theModificationTokens.size() )
//         {
//           return false;
//         }
      
//       for(std::vector<MolTokenStr>::const_iterator thisIter = theMolTokens.begin(), 
//             std::vector<MolTokenStr>::const_iterator otherIter = other.theMolTokens.begin();
//           thisIter != theMolTokens.end();
//           ++thisIter, ++otherIter)
//         {
//           if (*thisIter != *otherIter) return false;
//         }

//       for(std::vector<BindingTokenStr>::const_iterator thisIter = theBindingTokens.begin(), 
//             otherIter = other.theBindingTokens.begin();
//           thisIter != theBindingTokens.end();
//           ++thisIter, ++otherIter)
//         {
//           if (*thisIter != *otherIter) return false;
//         }

//       for(std::vector<ModificationTokenStr>::const_iterator thisIter = theModificationTokens.begin(), 
//             otherIter = other.theModificationTokens.begin();
//           thisIter != theModificationTokens.end();
//           ++thisIter, ++otherIter)
//         {
//           if (*thisIter != *otherIter) return false;
//         }
        
//       // If we've gone through all that and still haven't show them to be different, chances are good they're 
//       // the same...

//       return true;


            return (theMolTokens == other.theMolTokens &&
                    theBindingTokens == other.theBindingTokens &&
                    theModificationTokens == other.theModificationTokens);
                }
      
        void addMolTokenToOutputState(MolTokenStrCref aMolToken)
        {
            theMolTokens.push_back(aMolToken);
        }
        void addBindingTokenToOutputState(BindingTokenStrCref aBindingToken)
        {
            theBindingTokens.push_back(aBindingToken);
        }

        void addModificationTokenToOutputState(ModificationTokenStrCref aModificationToken)
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

std::ostream& operator<<(std::ostream& os, const nmr::ComplexOutputState& cos);
  
#endif
