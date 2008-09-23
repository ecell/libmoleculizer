/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2007, 2008 The Molecular Sciences Institute
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

#ifndef __ABSTRACTMOL_HH
#define __ABSTRACTMOL_HH

#include "nmrExceptions.hh"

#include "utl/macros.hh"

#include <string>
#include <vector>
#include <list>
#include <utility>

namespace nmr
{
    // An abstract class that implements a Molecule with a name.
    // Named molecules have names, binding sites, and modification sites, such that 
    // the binding sites can be either "bound" or "unbound", and modification sites
    // can be associated with any kind of modification.  

    DECLARE_CLASS( Mol );
    class Mol
    {
    public:

        DECLARE_TYPE(std::string, MolType);
        DECLARE_TYPE(std::string, BindingSite);
        DECLARE_TYPE(std::string, ModificationSite);
        DECLARE_TYPE(std::string, ModificationValue);

        DECLARE_TYPE(std::vector<BindingSite>, BindingSiteList);
        typedef std::pair<ModificationSite, ModificationValue> _ModificationToken;
        DECLARE_TYPE( _ModificationToken, ModificationToken);
        DECLARE_TYPE( std::vector<ModificationToken>,  ModificationList );

        Mol(MolTypeCref aMolType) 
            : 
            theMolType(aMolType) 
        {}
      
        virtual ~Mol()
        {}
      
        inline MolTypeCref 
        getMolType() const
        { 
            return theMolType; 
        }

        // API functions to build up the mol.
        virtual void addNewBindingSite( BindingSiteCref aBindingSite ) = 0;
        virtual void addNewModificationSite( ModificationSiteCref newModSite, ModificationValueCref modValue) = 0;

        virtual BindingSiteList getBindingList() const = 0;

        // To get the structure of the mol.
        virtual bool checkIfBindingSiteExists(BindingSiteCref aBindingSite) const = 0;
        virtual bool checkIfModificationSiteExists(ModificationSiteCref aModificationSite) const =0;
        virtual int getBindingSiteInteger(BindingSiteCref aBindingSite) const throw(NoSuchBindingSiteXcpt) =0; 
        virtual int getModificationSiteInteger(ModificationSiteCref aModificationSite) const throw(NoSuchModificationSiteXcpt) =0;
        virtual ModificationList getModificationList() const = 0;

        // To determine the state of the mol.
        virtual bool checkIfBindingSiteIsBound(BindingSiteCref aBindingSite) const throw(NoSuchBindingSiteXcpt) =0;
        virtual ModificationValue getModificationValueAtModificationSite(ModificationSiteCref aModificationSite) const throw(NoSuchModificationSiteXcpt) =0;        


        // To change the state of the mol.
        virtual void bindAtBindingSite(BindingSiteCref aBindingSite) throw(NoSuchBindingSiteXcpt) =0;
        virtual void unbindAtBindingSite(BindingSiteCref aBindingSite) throw(NoSuchBindingSiteXcpt) =0;
        virtual void updateModificationState(ModificationSiteCref aModificationSite, ModificationValueCref aModificationValue) throw(NoSuchModificationSiteXcpt) =0;

        virtual bool 
        operator<(MolCref rhsMol) const
        {
            return getMolType() < rhsMol.getMolType();
        }

    protected:
        const MolType theMolType;
    };
}

#endif 
