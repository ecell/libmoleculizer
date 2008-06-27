/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2008  Nathan Addy
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


#include <ostream>
#include "complexOutputState.hh"

namespace nmr
{
    bool 
    ComplexOutputState::operator==(ComplexOutputStateCref other) const
    {
        return (theMolTokens == other.theMolTokens &&
                theBindingTokens == other.theBindingTokens &&
                theModificationTokens == other.theModificationTokens);
    }

    bool 
    ComplexOutputState::operator!=(ComplexOutputStateCref other) const
    {
        return !( *this == other );
    }
      
    void 
    ComplexOutputState::addMolTokenToOutputState(MolTokenStrCref aMolToken)
    {
        theMolTokens.push_back(aMolToken);
    }

    void 
    ComplexOutputState::addBindingTokenToOutputState(BindingTokenStrCref aBindingToken)
    {
        theBindingTokens.push_back(aBindingToken);
    }

    void 
    ComplexOutputState::addModificationTokenToOutputState(ModificationTokenStrCref aModificationToken)
    {
        theModificationTokens.push_back(aModificationToken);
    }

    void 
    ComplexOutputState::clear()
    {
        theMolTokens.clear();
        theBindingTokens.clear();
        theModificationTokens.clear();
    }
}

std::ostream& operator<<(std::ostream& os, const nmr::ComplexOutputState& cos)
{
  os << "Mols:\n";
  for(std::vector<nmr::ComplexOutputState::MolTokenStr>::const_iterator iter = cos.theMolTokens.begin();
      iter != cos.theMolTokens.end();
      ++iter)
    {
      os << *iter << ", ";
          
    }
  os << "\nBindings:\n";

  for(std::vector<nmr::ComplexOutputState::BindingTokenStr>::const_iterator iter = cos.theBindingTokens.begin();
      iter != cos.theBindingTokens.end();
      ++iter)
    {
      os << "((" << (*iter).first.first << ", " << (*iter).first.second << "), (" << (*iter).second.first << ", " << (*iter).second.second << ")), ";
    }
      
  os << "\nModifications:\n";

  for(std::vector<nmr::ComplexOutputState::ModificationTokenStr>::const_iterator iter = cos.theModificationTokens.begin();
      iter != cos.theModificationTokens.end();
      ++iter)
    {
      os << "( " << (*iter).first << ", ( " << (*iter).second.first << ", " << (*iter).second.second << ")), ";
    }
  if (cos.theModificationTokens.size() == 0) os << "*";

  os << "\n";
  return os;
}
