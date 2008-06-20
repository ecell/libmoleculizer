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


#ifndef __MINIMALMOL_HH
#define __MINIMALMOL_HH

#include "abstractMol.hh"

#include <map>
#include <string>
#include <vector>
#include <set>

namespace nmr
{

    // This class is a "MinimalMol" because it has no a priori structure.
    // BindingSites and ModificationSites are always presumed to exist.

    DECLARE_CLASS( MinimalMol );
    class MinimalMol 
        : 
        public Mol
    {
    public:

        MinimalMol( MolTypeCref molType)
            :
            Mol( molType)
        {}

        MinimalMol( MolTypeCref molType, BindingSiteListCref bindingSitesList ) 
            : 
            Mol( molType) 
        {
//             for(BindingSiteList::const_iterator iter = bindingSitesList.begin();
//                 iter != bindingSitesList.end();
//                 ++iter)
//             {
//                 theBindingSiteStates.insert( std::make_pair(*iter,
//                                                             false) );
//             }
        }

        bool checkIfBindingSiteExists(BindingSiteCref aBindingSite) const;
        bool checkIfModificationSiteExists( ModificationSiteCref aModificationSite) const;

        bool 
        checkIfBindingSiteIsBound(BindingSiteCref aBindingSite) 
            const throw(NoSuchBindingSiteXcpt);

        void 
        addNewBindingSite( BindingSiteCref aBindingSite);


        ModificationValue 
        getModificationValueAtModificationSite(ModificationSiteCref aModificationSite) 
            const throw(NoSuchModificationSiteXcpt);

        ModificationList 
        getModificationList() const;

        void 
        bindAtBindingSite(BindingSiteCref aBindingSite) 
            throw(NoSuchBindingSiteXcpt);

        void 
        unbindAtBindingSite(BindingSiteCref aBindingSite) 
            throw(NoSuchBindingSiteXcpt);

        void 
        updateModificationState(ModificationSiteCref aModificationSite,
                                ModificationValueCref aModificationValue) 
            throw(NoSuchModificationSiteXcpt);
    
        int 
        getBindingSiteInteger(BindingSiteCref aBindingSite) 
            const throw(NoSuchBindingSiteXcpt);

        int 
        getModificationSiteInteger(ModificationSiteCref aModificationSite) 
            const throw(NoSuchModificationSiteXcpt);    

        void
        addNewModificationSite( ModificationSiteCref newModSite,
                                ModificationValueCref modValue);
//    protected:
        // std::map<BindingSite, bool> theBindingSiteStates; // true means bound, false,unbound.

        std::map<BindingSite, unsigned int> bindingSiteNameToNdxMap;

        // This is maybe not the best... Dragons! Beware!
        std::vector<bool> bindingSiteIsBound;

        std::map<ModificationSite, ModificationValue> theModificationStates;
        std::map<ModificationSite, std::set<ModificationValue> > theLegalModifications;
    };

}

#endif

