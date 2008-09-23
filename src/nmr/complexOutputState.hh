/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2007, 2008 The Molecular Sciences Institute
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
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
        DECLARE_TYPE( std::string, MolTokenStr);
        typedef std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > __BindingTokenStr;
        typedef std::pair<std::string, std::pair<std::string, std::string> > __ModificationTokenStr;

        DECLARE_TYPE( __BindingTokenStr, BindingTokenStr);
        DECLARE_TYPE( __ModificationTokenStr, ModificationTokenStr);

        std::vector<MolTokenStr> theMolTokens;
        std::vector<BindingTokenStr> theBindingTokens;
        std::vector<ModificationTokenStr> theModificationTokens;

        bool operator==(const ComplexOutputState& other) const;
        bool operator!=(const ComplexOutputState& other) const;
        
        void addMolTokenToOutputState(MolTokenStrCref aMolToken);
        void addBindingTokenToOutputState(BindingTokenStrCref aBindingToken);
        void addModificationTokenToOutputState(ModificationTokenStrCref aModificationToken);

        void clear();
        std::string repr() const;
    };
}

std::ostream& operator<<(std::ostream& os, const nmr::ComplexOutputState& cos);
  
#endif
