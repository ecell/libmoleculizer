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

#include "nmr/namedMolecule.hh"

namespace nmr
{
    void 
    MinimalMol::addNewBindingSite( BindingSiteCref aBindingSite)
    {
        if (theBindingSiteStates.find(aBindingSite) != theBindingSiteStates.end())
        {
            throw GeneralNmrXcpt( "Binding Site already has been added in MinimalMol::addNewBindingSite. (key: lfjdkas)");
        }

        theBindingSiteStates.insert( std::make_pair( aBindingSite, false) );
    }

    bool 
    MinimalMol::checkIfBindingSiteExists( BindingSiteCref aBindingSite) const
    {
        return (theBindingSiteStates.find(aBindingSite) != theBindingSiteStates.end() );
    }

    bool
    MinimalMol::checkIfModificationSiteExists( ModificationSiteCref aModificationSite) const
    {
        return (theModificationStates.find( aModificationSite) != theModificationStates.end() );
    }

    bool 
    MinimalMol::checkIfBindingSiteIsBound(BindingSiteCref aBindingSite) 
        const throw(NoSuchBindingSiteXcpt)
    {
        std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            throw NoSuchBindingSiteXcpt( getMolType(),
                                         aBindingSite);
        }
        else
        {
            return aBindingSiteLocation->second;
        }
    }

    MinimalMol::ModificationValue 
    MinimalMol::getModificationValueAtModificationSite(ModificationSiteCref aModificationSite) 
        const throw(NoSuchModificationSiteXcpt)
    {
        std::map<ModificationSite, ModificationValue>::const_iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
        if (aModificationSiteLocation == theModificationStates.end())
        {
            // This apparently violates the T&C of MinimalMol, but I can't think of anyway to 
            // handle this.  Of course, any subsequent call to 
            // updateModificationValueAtModificationSite(aModificationSite, foo) would be successful.
            throw NoSuchModificationSiteXcpt( getMolType(), aModificationSite );
        }
        
        return aModificationSiteLocation->second;
    }

    MinimalMol::ModificationList 
    MinimalMol::getModificationList() const
    {
        return ModificationList( theModificationStates.begin(), theModificationStates.end());
    }

    void 
    MinimalMol::bindAtBindingSite(BindingSiteCref aBindingSite) 
        throw(NoSuchBindingSiteXcpt)
    {

        if (theBindingSiteStates.find(aBindingSite) == theBindingSiteStates.end()) 
        {
            throw NoSuchBindingSiteXcpt( getMolType(), aBindingSite);
        }
        else if (theBindingSiteStates[aBindingSite] == true) 
        {
            // The only way this can fail is if the bindingSite is already bound.
            throw BindingSiteAlreadyBoundXcpt( getMolType(), aBindingSite );
        }
        else
        {
            theBindingSiteStates[aBindingSite] = true;
        }

    }

    void 
    MinimalMol::unbindAtBindingSite(BindingSiteCref aBindingSite) 
        throw(NoSuchBindingSiteXcpt)
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);

        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            throw NoSuchBindingSiteXcpt( getMolType(), aBindingSite);
            
        }
        else if (theBindingSiteStates[aBindingSite] == false)
        {
            throw BindingSiteAlreadyUnboundXcpt( getMolType(), aBindingSite);
        }
        else
        {
            theBindingSiteStates[aBindingSite] = false;
        }
   
    }

    void
    MinimalMol::updateModificationState(ModificationSiteCref aModificationSite,
                                        ModificationValueCref aModificationValue)
        throw(NoSuchModificationSiteXcpt)
    {
        std::map<ModificationSite, ModificationValue>::iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
        if (aModificationSiteLocation == theModificationStates.end())
        {
            theModificationStates.insert( std::make_pair( aModificationSite,
                                                          aModificationValue ));
        }
        else
        {
            theModificationStates[ aModificationSite ] = aModificationValue;
        }

    }

    int 
    MinimalMol::getBindingSiteInteger(BindingSiteCref aBindingSite) 
        const throw(NoSuchBindingSiteXcpt)
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            throw NoSuchBindingSiteXcpt( getMolType(), aBindingSite);
        }

        int index = 0;
        for(std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.begin();
            aBindingSiteLocation != theBindingSiteStates.end();
            ++aBindingSiteLocation, ++index)
        {
            if (aBindingSiteLocation->first == aBindingSite)
            {
                return index;
            }
        }
        
        throw utl::xcpt( "Unexpected Error in MinimalMol::getBindingSiteInteger. Please contact the maintainer of this code.");
    }

    int 
    MinimalMol::getModificationSiteInteger(ModificationSiteCref aModificationSite) 
        const throw(NoSuchModificationSiteXcpt)
    {
        std::map<ModificationSite, ModificationValue>::const_iterator aModSiteLoc = theModificationStates.find(aModificationSite);
        if (aModSiteLoc == theModificationStates.end())
        {
            throw NoSuchModificationSiteXcpt( getMolType(), aModificationSite);
        }
    
        int index = 0;
        for(std::map<ModificationSite, ModificationValue>::const_iterator aModificationSiteLocationIter = theModificationStates.begin();
            aModificationSiteLocationIter != theModificationStates.end();
            ++aModificationSiteLocationIter, ++index)
        {
            if (aModificationSiteLocationIter->first == aModificationSite) return index;

        }

        throw utl::xcpt( "Unexpected Error in MinimalMol::getModificationSiteInteger. Please contact the maintainer of this code.");
    }


    void
    MinimalMol::addNewModificationSite( ModificationSiteCref newModSite,
                                        ModificationValueCref modValue)
    {
        if (theModificationStates.find( newModSite ) != theModificationStates.end() )
        {
            throw GeneralNmrXcpt( "Error, modification site already exists and so cannot be added anew in MinimalMol::addNewModificationSite.  (key: kjdfha)");
        }
        
        theModificationStates.insert( std::make_pair( newModSite,
                                                      modValue ) );
    }
    


}

