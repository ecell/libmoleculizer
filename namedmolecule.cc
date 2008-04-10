#include "nmr/mol.hh"

namespace nmr
{

    bool 
    Mol::checkIfBindingSiteExists(BindingSite aBindingSite) const 
    {
        std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        return ( aBindingSiteLocation != theBindingSiteStates.end() );
    }

    bool 
    Mol::checkIfModificationSiteExists(ModificationSite aModificationSite) const
    {
        std::map<ModificationSite, ModificationValue>::const_iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
        return ( aModificationSiteLocation != theModificationStates.end() );
    }

    bool 
    Mol::checkIfBindingSiteIsBound(BindingSite aBindingSite) const
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            throw CSXcpt("Mol::checkIfBindingSiteIsBound", "There is no BindingSite named " + aBindingSite+".");
        }
 
        return aBindingSiteLocation->second;
    }

  
    Mol::ModificationValue 
    Mol::getModificationValueAtModificationSite(ModificationSite aModificationSite) const
    {
        //Make sure the ModificationSite exists.
        std::map<ModificationSite, ModificationValue>::const_iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
        if (aModificationSiteLocation == theModificationStates.end())
        {
            throw CSXcpt("Mol::getModificationValueAtModificationSite", "There is no such ModificationSite as "+aModificationSite);
        }

        return aModificationSiteLocation->second;
    }

    Mol::ModificationList 
    Mol::getModificationList() const
    {
        return ModificationList( theModificationStates.begin(), theModificationStates.end());
    
    }


    void 
    Mol::bindAtBindingSite(BindingSite aBindingSite)
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            throw CSXcpt("Mol::bindAtBindingSite", "There is no BindingSite named " + aBindingSite+".");
        }

        if (aBindingSiteLocation->second == true)
        {
            std::string anErrorMessage = "Error: Mol::bindAtBindingSite failed.\n"+aBindingSite+"is already bound.";
            throw CSXcpt(anErrorMessage);
        }
    
        theBindingSiteStates[aBindingSite] = true;
   
    }


    void 
    Mol::unbindAtBindingSite(BindingSite aBindingSite)
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            std::string anErrorMessage = "Error: Mol::unbindAtBindingSite failed.\nThere is no BindingSite named " + aBindingSite+".";
            throw CSXcpt(anErrorMessage);
        }

        if (aBindingSiteLocation->second == false)
        {
            std::string anErrorMessage = "Error: Mol::unbindAtBindingSite failed.\n"+aBindingSite+"is already unbound.";
            throw CSXcpt(anErrorMessage);
        }
    
        theBindingSiteStates[aBindingSite] = false;
   
    }

    void
    Mol::updateModificationState(ModificationSite aModificationSite,
                                 ModificationValue aModificationValue)
    {
        //Make sure the ModifcationSite exists
        std::map<ModificationSite, ModificationValue>::iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
        if (aModificationSiteLocation == theModificationStates.end())
        {
            std::string anErrorMessage = "Error: Mol::updateModification failed.\nThere is no such ModificationSite as "+aModificationSite;
            throw CSXcpt(anErrorMessage);
        }
    
        //Make sure the ModificationValue is legal for this site.
        std::set<ModificationValue>::iterator aModificationValueLocation = theLegalModifications[aModificationSite].find(aModificationValue);
        if ( aModificationValueLocation == theLegalModifications[aModificationSite].end() )
        {
            std::string anErrorMessage = "Error: Mol::updateModificationState failed.\n" + aModificationValue + " is not a legal ModificationValue at this ModifcationSite";
            throw CSXcpt(anErrorMessage);
        }
    
        theModificationStates[ aModificationSite ] = aModificationValue;
    }



    int 
    Mol::getBindingSiteInteger(BindingSite aBindingSite) const
    {
        // Make sure the site exists.
        std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
        if (aBindingSiteLocation == theBindingSiteStates.end()) 
        {
            throw CSXcpt("Mol::bindAtBindingSite", "There is no BindingSite named " + aBindingSite+".");
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
        throw CSXcpt("Mol::getBindingSiteInteger(BindingSite aBindingSite) const", "aBindingSite should have been found by iterating through theBindingSiteStates but wasn't.");
    }

    int 
    Mol::getModificationSiteInteger(ModificationSite aModificationSite) const
    {
        std::map<ModificationSite, ModificationValue>::const_iterator aModSiteLoc = theModificationStates.find(aModificationSite);
        if (aModSiteLoc == theModificationStates.end())
        {
            throw CSXcpt("Mol::getModificationSiteInteger(ModificationSite aModificationSite) const", "aModificationSite was not found in theModificationStates.");
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
        throw CSXcpt("Mol::getModificationSiteInteger(ModificationSite aModificationSite) const", "aModificationSite was not found by iterating through theModificationState");
    }

    void 
    Mol::addNewBindingSite(BindingSite aBindingSite)
    {

        // Check if the Binding site already exists.
        // If it does, throw and exception
        std::map<BindingSite, bool>::iterator loc = theBindingSiteStates.find(aBindingSite);
        if (loc != theBindingSiteStates.end())
        {
            std::string anErrorMessage = "Error: Mol::addNewBindingSite failed.\nBindingSite "+aBindingSite+" already occurs in the list of BindingSites.";
            throw CSXcpt(anErrorMessage);
        }

        // Add the BindingSite into the map of binding sites with an unbound state.
        theBindingSiteStates.insert( std::make_pair(aBindingSite, false) );
    }

    void 
    Mol::addNewModificationSite(ModificationSite aModificationSite,
                                ListOfModificationValues aListOfValidModificationValues)
    {
        // Make sure the Modification site does not already exist.
        std::map<ModificationSite, ModificationValue>::iterator loc = theModificationStates.find(aModificationSite);
        if (loc != theModificationStates.end())
        {
            std::string anErrorMessage = "Error: Mol::addNewModificationSite failed.\nModicationSite "+aModificationSite+" already occurs in the list of ModificationSites.";
            throw CSXcpt(anErrorMessage);
        }

        if(aListOfValidModificationValues.empty())
        {
            throw CSXcpt("Error: Mol::addNewModificationSite failed.\n listOfValidModficationTypes was empty.");
        }

        // Check if the listOfValidModifications are all unique.  If not, throw an exception.
        // This is strictly not necessary, but "better safe than sorry" seems to apply here.
        std::set<ModificationValue> uniqueModificationValues(aListOfValidModificationValues.begin(), aListOfValidModificationValues.end());
        if (uniqueModificationValues.size() != aListOfValidModificationValues.size())
        {
            std::string anErrorMessage = "Error: Mol::addNewModificationSite failed.\nThe list of ModificationValues were not all unique.";
            throw CSXcpt(anErrorMessage);
        }

        // Add the ModificationSite along with the full set of ModificationValues into theLegalModifications.
        theLegalModifications.insert( std::make_pair(aModificationSite, uniqueModificationValues));
    
        //Add the ModificationSite along with its default modification into theDefaultModifications.
        theDefaultModificationValues.insert( std::make_pair(aModificationSite, aListOfValidModificationValues[0]) );

        //Add the ModificationSite along with its default value into theModificationStates.
        theModificationStates.insert( std::make_pair(aModificationSite, theDefaultModificationValues[aModificationSite]) );
    }

    void 
    Mol::addNewModificationSite(ModificationSite aModificationSite)
    {
        // Make sure the Modification site does not already exist.
        std::map<ModificationSite, ModificationValue>::iterator loc = theModificationStates.find(aModificationSite);
        if (loc != theModificationStates.end())
        {
            std::string anErrorMessage = "Error: Mol::addNewModificationSite failed.\nModicationSite "+aModificationSite+" already occurs in the list of ModificationSites.";
            throw CSXcpt(anErrorMessage);
        }

        // Add the ModificationSite.
        theLegalModifications.insert( std::make_pair(aModificationSite, uniqueModificationValues));






    }

}

