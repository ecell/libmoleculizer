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

#include "nmr/namedmolecule.hh"

namespace nmr
{
  
  bool 
  MinimalMol::checkIfBindingSiteExists(BindingSite aBindingSite) 
    const
  {
    std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
    return ( aBindingSiteLocation != theBindingSiteStates.end() );
  }

  bool 
  MinimalMol::checkIfModificationSiteExists(ModificationSite aModificationSite) 
    const 
  {
    std::map<ModificationSite, ModificationValue>::const_iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
    return ( aModificationSiteLocation != theModificationStates.end() );
  }

  bool 
  MinimalMol::checkIfBindingSiteIsBound(BindingSite aBindingSite) 
    const throw(NoSuchBindingSiteXcpt)
  {
    // Make sure the site exists.
    std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
    if (aBindingSiteLocation == theBindingSiteStates.end()) 
      {
        throw CSXcpt("Mol::checkIfBindingSiteIsBound", "There is no BindingSite named " + aBindingSite+".");
      }
 
    return aBindingSiteLocation->second;
  }

  MinimalMol::ModificationValue 
  MinimalMol::getModificationValueAtModificationSite(ModificationSite aModificationSite) 
    const throw(NoSuchModificationSiteXcpt)
  {
    //Make sure the ModificationSite exists.
    std::map<ModificationSite, ModificationValue>::const_iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
    if (aModificationSiteLocation == theModificationStates.end())
      {
        throw CSXcpt("MinimalMol::getModificationValueAtModificationSite", "There is no such ModificationSite as "+aModificationSite);
      }

    return aModificationSiteLocation->second;
  }

  MinimalMol::ModificationList 
  MinimalMol::getModificationList() const
  {
    return ModificationList( theModificationStates.begin(), theModificationStates.end());
    
  }


  void 
  MinimalMol::bindAtBindingSite(BindingSite aBindingSite) 
    throw(NoSuchBindingSiteXcpt)
  {
    // Make sure the site exists.
    std::map<BindingSite, bool>::iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
    if (aBindingSiteLocation == theBindingSiteStates.end()) 
      {
        throw CSXcpt("MinimalMol::bindAtBindingSite", "There is no BindingSite named " + aBindingSite+".");
      }

    if (aBindingSiteLocation->second == true)
      {
        std::string anErrorMessage = "Error: MinimalMol::bindAtBindingSite failed.\n"+aBindingSite+"is already bound.";
        throw CSXcpt(anErrorMessage);
      }
    
    theBindingSiteStates[aBindingSite] = true;
   
  }


  void 
  MinimalMol::unbindAtBindingSite(BindingSite aBindingSite) 
    throw(NoSuchBindingSiteXcpt)
  {
    // Make sure the site exists.
    std::map<BindingSite, bool>::iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
    if (aBindingSiteLocation == theBindingSiteStates.end()) 
      {
        std::string anErrorMessage = "Error: MinimalMol::unbindAtBindingSite failed.\nThere is no BindingSite named " + aBindingSite+".";
        throw CSXcpt(anErrorMessage);
      }

    if (aBindingSiteLocation->second == false)
      {
        std::string anErrorMessage = "Error: MinimalMol::unbindAtBindingSite failed.\n"+aBindingSite+"is already unbound.";
        throw CSXcpt(anErrorMessage);
      }
    
    theBindingSiteStates[aBindingSite] = false;
   
  }

  void
  MinimalMol::updateModificationState(ModificationSite aModificationSite,
                                      ModificationValue aModificationValue)
    throw(NoSuchModificationSiteXcpt)
  {
    //Make sure the ModifcationSite exists
    std::map<ModificationSite, ModificationValue>::iterator aModificationSiteLocation = theModificationStates.find(aModificationSite);
    if (aModificationSiteLocation == theModificationStates.end())
      {
        std::string anErrorMessage = "Error: MinimalMol::updateModification failed.\nThere is no such ModificationSite as "+aModificationSite;
        throw CSXcpt(anErrorMessage);
      }
    
    //Make sure the ModificationValue is legal for this site.
    std::set<ModificationValue>::iterator aModificationValueLocation = theLegalModifications[aModificationSite].find(aModificationValue);
    if ( aModificationValueLocation == theLegalModifications[aModificationSite].end() )
      {
        std::string anErrorMessage = "Error: MinimalMol::updateModificationState failed.\n" + aModificationValue + " is not a legal ModificationValue at this ModifcationSite";
        throw CSXcpt(anErrorMessage);
      }
    
    theModificationStates[ aModificationSite ] = aModificationValue;
  }



  int 
  MinimalMol::getBindingSiteInteger(BindingSite aBindingSite) 
    const throw(NoSuchBindingSiteXcpt)
  {
    // Make sure the site exists.
    std::map<BindingSite, bool>::const_iterator aBindingSiteLocation = theBindingSiteStates.find(aBindingSite);
    if (aBindingSiteLocation == theBindingSiteStates.end()) 
      {
        throw CSXcpt("MinimalMol::bindAtBindingSite", "There is no BindingSite named " + aBindingSite+".");
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
    throw CSXcpt("MinimalMol::getBindingSiteInteger(BindingSite aBindingSite) const", "aBindingSite should have been found by iterating through theBindingSiteStates but wasn't.");
  }

  int 
  MinimalMol::getModificationSiteInteger(ModificationSite aModificationSite) 
    const throw(NoSuchModificationSiteXcpt)
  {
    std::map<ModificationSite, ModificationValue>::const_iterator aModSiteLoc = theModificationStates.find(aModificationSite);
    if (aModSiteLoc == theModificationStates.end())
      {
        throw CSXcpt("MinimalMol::getModificationSiteInteger(ModificationSite aModificationSite) const", "aModificationSite was not found in theModificationStates.");
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
    throw CSXcpt("MinimalMol::getModificationSiteInteger(ModificationSite aModificationSite) const", "aModificationSite was not found by iterating through theModificationState");
  }

  void 
  MinimalMol::addNewBindingSite(BindingSite aBindingSite)
    throw(NoSuchBindingSiteXcpt)
  {

    // Check if the Binding site already exists.
    // If it does, throw and exception
    std::map<BindingSite, bool>::iterator loc = theBindingSiteStates.find(aBindingSite);
    if (loc != theBindingSiteStates.end())
      {
        std::string anErrorMessage = "Error: MinimalMol::addNewBindingSite failed.\nBindingSite "+aBindingSite+" already occurs in the list of BindingSites.";
        throw CSXcpt(anErrorMessage);
      }

    // Add the BindingSite into the map of binding sites with an unbound state.
    theBindingSiteStates.insert( std::make_pair(aBindingSite, false) );
  }

  void 
  MinimalMol::addNewModificationSite(ModificationSite aModificationSite,
                                     ListOfModificationValues aListOfValidModificationValues)
    throw(NoSuchModificationSiteXcpt)
  {
    // Make sure the Modification site does not already exist.
    std::map<ModificationSite, ModificationValue>::iterator loc = theModificationStates.find(aModificationSite);
    if (loc != theModificationStates.end())
      {
        std::string anErrorMessage = "Error: MinimalMol::addNewModificationSite failed.\nModicationSite "+aModificationSite+" already occurs in the list of ModificationSites.";
        throw CSXcpt(anErrorMessage);
      }

    if(aListOfValidModificationValues.empty())
      {
        throw CSXcpt("Error: MinimalMol::addNewModificationSite failed.\n listOfValidModficationTypes was empty.");
      }

    // Check if the listOfValidModifications are all unique.  If not, throw an exception.
    // This is strictly not necessary, but "better safe than sorry" seems to apply here.
    std::set<ModificationValue> uniqueModificationValues(aListOfValidModificationValues.begin(), aListOfValidModificationValues.end());
    if (uniqueModificationValues.size() != aListOfValidModificationValues.size())
      {
        std::string anErrorMessage = "Error: MinimalMol::addNewModificationSite failed.\nThe list of ModificationValues were not all unique.";
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
  MinimalMol::addNewModificationSite(ModificationSite aModificationSite)
    throw(NoSuchModificationSiteXcpt)
  {
    // Make sure the Modification site does not already exist.
    std::map<ModificationSite, ModificationValue>::iterator loc = theModificationStates.find(aModificationSite);
    if (loc != theModificationStates.end())
      {
        std::string anErrorMessage = "Error: MinimalMol::addNewModificationSite failed.\nModicationSite "+aModificationSite+" already occurs in the list of ModificationSites.";
        throw CSXcpt(anErrorMessage);
      }

    // Add the ModificationSite.
    theLegalModifications.insert( std::make_pair(aModificationSite, uniqueModificationValues));
  }

}

