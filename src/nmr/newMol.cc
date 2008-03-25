#include "nmr/newMol.hh"

namespace nmr
{

    bool 
    SimpleMol::checkIfBindingSiteExists(BindingSite aBindingSite) const 
    {
        std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        return ( aBindingSiteLocation != theBindingSiteStates.end() );
    }

    bool 
    SimpleMol::checkIfModificationSiteExists(ModificationSite aModificationSite) const
    {
        std::map<ModificationSite, ModificationValue>::const_iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
        return ( aModificationSiteLocation != theModificationStates.end() );
    }

    bool 
    SimpleMol::checkIfBindingSiteIsBound(BindingSite aBindingSite) const
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            throw CSXcpt("SimpleMol::checkIfBindingSiteIsBound", "There is no BindingSite named " + aBindingSite+".");
        }
 
        return aBindingSiteLocation->second;
    }

  
    SimpleMol::ModificationValue 
    SimpleMol::getModificationValueAtModificationSite(ModificationSite aModificationSite) const
    {
        //Make sure the ModificationSite exists.
        std::map<ModificationSite, ModificationValue>::const_iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
        if (aModificationSiteLocation == theModificationStates.end())
        {
            throw CSXcpt("SimpleMol::getModificationValueAtModificationSite", "There is no such ModificationSite as "+aModificationSite);
        }

        return aModificationSiteLocation->second;
    }

    SimpleMol::ModificationList 
    SimpleMol::getModificationList() const
    {
        return ModificationList( theModificationStates.begin(), 
                                 theModificationStates.end() );
    
    }


    void 
    SimpleMol::bindAtBindingSite(BindingSite aBindingSite)
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            throw CSXcpt("SimpleMol::bindAtBindingSite", "There is no BindingSite named " + aBindingSite+".");
        }

        if (aBindingSiteLocation->second == true)
        {
            std::string anErrorMessage = "Error: SimpleMol::bindAtBindingSite failed.\n"+aBindingSite+"is already bound.";
            throw CSXcpt(anErrorMessage);
        }
    
        theBindingSiteStates[aBindingSite] = true;
   
    }


    void 
    SimpleMol::unbindAtBindingSite(BindingSite aBindingSite)
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            std::string anErrorMessage = "Error: SimpleMol::unbindAtBindingSite failed.\nThere is no BindingSite named " + aBindingSite+".";
            throw CSXcpt(anErrorMessage);
        }

        if (aBindingSiteLocation->second == false)
        {
            std::string anErrorMessage = "Error: SimpleMol::unbindAtBindingSite failed.\n"+aBindingSite+"is already unbound.";
            throw CSXcpt(anErrorMessage);
        }
    
        theBindingSiteStates[aBindingSite] = false;
   
    }

    void
    SimpleMol::updateModificationState(ModificationSite aModificationSite,
                                       ModificationValue aModificationValue)
    {
        //Make sure the ModifcationSite exists
        std::map<ModificationSite, ModificationValue>::iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
        if (aModificationSiteLocation == theModificationStates.end())
        {
            std::string anErrorMessage = "Error: SimpleMol::updateModification failed.\nThere is no such ModificationSite as "+aModificationSite;
            throw CSXcpt(anErrorMessage);
        }
    
        theModificationStates[ aModificationSite ] = aModificationValue;
    }



    int 
    SimpleMol::getBindingSiteInteger(BindingSite aBindingSite) const
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            throw CSXcpt("SimpleMol::bindAtBindingSite", "There is no BindingSite named " + aBindingSite+".");
        }

        int index=0;
        for(std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.begin();
            aBindingSiteLocation != theBindingSiteStates.end();
            ++aBindingSiteLocation, ++index)
        {
            if (aBindingSiteLocation->first == aBindingSite)
            {
                return index;
            }
        }
    
        // Error, we should have picked it up.
        throw CSXcpt("SimpleMol::getBindingSiteInteger(BindingSite aBindingSite) const", "aBindingSite should have been found by iterating through theBindingSiteStates but wasn't.");
    }

    int 
    SimpleMol::getModificationSiteInteger(ModificationSite aModificationSite) const
    {
        std::map<ModificationSite, ModificationValue>::const_iterator aModSiteLoc = theModificationStates.find(aModificationSite);
        if (aModSiteLoc == theModificationStates.end())
        {
            throw CSXcpt("SimpleMol::getModificationSiteInteger(ModificationSite aModificationSite) const", "aModificationSite was not found in theModificationStates.");
        }
    
        int index = 0;
        for(std::map<ModificationSite, ModificationValue>::const_iterator aModificationSiteLocationIter = theModificationStates.begin();
            aModificationSiteLocationIter != theModificationStates.end();
            ++aModificationSiteLocationIter, ++index)
        {
            if (aModificationSiteLocationIter->first == aModificationSite)
            {
                return index;
            }
        }
        throw CSXcpt("SimpleMol::getModificationSiteInteger(ModificationSite aModificationSite) const", "aModificationSite was not found by iterating through theModificationState");
    }

    void 
    SimpleMol::addNewBindingSite(BindingSite aBindingSite)
    {

        // Check if the Binding site already exists.
        // If it does, throw and exception
        std::map<BindingSite, bool>::iterator loc = theBindingSiteStates.find(aBindingSite);
        if (loc != theBindingSiteStates.end())
        {
            std::string anErrorMessage = "Error: SimpleMol::addNewBindingSite failed.\nBindingSite "+aBindingSite+" already occurs in the list of BindingSites.";
            throw CSXcpt(anErrorMessage);
        }

        // Add the BindingSite into the map of binding sites with an unbound state.
        theBindingSiteStates.insert( std::make_pair(aBindingSite, false) );
    }

    void 
    SimpleMol::addNewModificationSite(ModificationSite aModificationSite, ModificationValue theModificationValue)
    {
        // Make sure the Modification site does not already exist.
        std::map<ModificationSite, ModificationValue>::iterator loc = theModificationStates.find(aModificationSite);
        if (loc != theModificationStates.end())
        {
            std::string anErrorMessage = "Error: SimpleMol::addNewModificationSite failed.\nModicationSite "+aModificationSite+" already occurs in the list of ModificationSites.";
            throw CSXcpt(anErrorMessage);
        }

        //Add the ModificationSite along with its default value into theModificationStates.
        theModificationStates.insert( std::make_pair(aModificationSite, theModificationValue) );
    }

}

