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


#ifndef NMR_MOL_HH
#define NMR_MOL_HH

#include "abstractMol.hh"
#include <map>
#include <string>
#include <vector>
#include <set>

namespace nmr
{
  class SimpleMol 
    : 
    public Mol
  {
    
  public:
    typedef std::vector<std::pair<ModificationSite, ModificationValue>  > ModificationList;
    
    SimpleMol( MolType theMolType ) 
      : 
      AbstractMol(theMolType)
    {
      ; // do nothing
    }
    
    bool checkIfBindingSiteExists(BindingSite aBindingSite) const;
    bool checkIfModificationSiteExists(ModificationSite aModificationSite) const;

    bool checkIfBindingSiteIsBound(BindingSite aBindingSite) const;

    ModificationValue getModificationValueAtModificationSite(ModificationSite aModificationSite) const;

    void bindAtBindingSite(BindingSite aBindingSite);
    void unbindAtBindingSite(BindingSite aBindingSite);
    void updateModificationState(ModificationSite aModificationSite,
                                 ModificationValue aModificationValue);
    
    ModificationList getModificationList() const;
    int getBindingSiteInteger(BindingSite aBindingSite) const;
    int getModificationSiteInteger(ModificationSite aModificationSite) const;

    void addNewBindingSite(BindingSite aBindingSite);
    void addNewModificationSite(ModificationSite , ModificationValue);

  protected:
    std::map<BindingSite, bool> theBindingSiteStates; // true means bound, false,unbound.
    std::map<ModificationSite, ModificationValue> theModificationStates;

  };


}

#endif
