/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
**
**        This file is part of Libmoleculizer
**
**        Copyright (C) 2001-2008 The Molecular Sciences Institute.
**
**::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
**
** Moleculizer is free software; you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published
** by the Free Software Foundation; either version 3 of the License, or
** (at your option) any later version.
**
** Moleculizer is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with Moleculizer; if not, write to the Free Software Foundation
** Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
**
** END HEADER
**
** Original Author:
**   Nathan Addy, Scientific Programmer, Molecular Sciences Institute, 2008
**
** Modifing Authors:
**
*/

#include "libmzr_c_interface.h"

#include "utl/xcpt.hh"
#include "moleculizer.hh"

// 
// Declarations for functions locally defined.
// 
// 

mzr::moleculizer* 
convertCMzrPtrToMzrPtr(moleculizer* cMzrPtr);

int 
getSpeciesMapSize( const mzr::mzrReaction::multMap& speciesMap);

species* 
createNewCSpeciesFromMzrSpecies( moleculizer* cMzrPtr, const mzr::mzrSpecies* pMzrSpecies);

reaction* 
createNewCRxnFromMzrReaction( moleculizer* cMzrPtr, const mzr::mzrReaction* pMzrReaction);

int 
convertUserNameToSpeciesKey( moleculizer* handle, char* theUserName, char* correspondingSpeciesKey);
 

void 
constructCSpeciesArrayFromSpeciesMap(moleculizer* cMzrPtr, 
                                     const mzr::mzrReaction::multMap& speciesMap, 
                                     species* speciesList, int& numberInList);



//
// The interface presented in the header file -- the c-interface.
//
// 


moleculizer* createNewMoleculizerObject()
{
// Creates a new moleculizer object (the c-struct) that wraps a mzr::moleculizer c++ object.
// The pointer is owned by whoever calls this function and must be free'd using the 
// freeMoleculizerObject function.

    moleculizer_handle* newHandle = new moleculizer_handle;
    newHandle->mzrObject = (void*) new mzr::moleculizer;
    return newHandle;
}

int setupSpatialExtrapolation(moleculizer* handle)
{
    enum LOCAL_ERROR_TYPE { SUCCESS = 0,
                            UNKNOWN_ERROR = 1,
                            MODEL_ALREADY_LOADED = 2 };


    mzr::moleculizer* underlyingMoleculizerObject = convertCMzrPtrToMzrPtr(handle);

    if ( underlyingMoleculizerObject->getModelHasBeenLoaded() ) return MODEL_ALREADY_LOADED;

    try
    {
        underlyingMoleculizerObject->setRateExtrapolation( false );
//        underlyingMoleculizerObject->enableSpatialReactionNetworkGeneration();
        return SUCCESS;
    }
    catch(...)
    {
        return UNKNOWN_ERROR;
    }

}

int setupNonSpatialExtrapolation(moleculizer* handle, int useMassBasedRateExtrapolation)
{

    enum LOCAL_ERROR_TYPE { SUCCESS = 0,
                            UNKNOWN_ERROR = 1,
                            MODEL_ALREADY_LOADED = 2 };

    mzr::moleculizer* underlyingMoleculizerObject = convertCMzrPtrToMzrPtr(handle);
    
    if ( underlyingMoleculizerObject->getModelHasBeenLoaded() ) return MODEL_ALREADY_LOADED;

    try
    {
        underlyingMoleculizerObject->setRateExtrapolation( useMassBasedRateExtrapolation );
//        underlyingMoleculizerObject->enableSpatialReactionNetworkGeneration();
        return SUCCESS;
    }
    catch(...)
    {
        return UNKNOWN_ERROR;
    }
}

void freeMoleculizerObject( moleculizer* handle)
{
// This function frees the c-style moleculizer handle.
    delete convertCMzrPtrToMzrPtr(handle);
    delete handle;
}


int attachRulesFileToMoleculizerObject(moleculizer* handle, char* fileName)
{
// This function takes a string to a file containing an xml rules 
// description for the system and loads it into mzr::moleculizer.
// See documentation for a specification of model files as well as
// examples.

   enum LOCAL_ERROR_TYPE { SUCCESS = 0,
                           UNKNOWN_ERROR = 1, 
                           DOCUMENT_UNPARSABLE = 2,
                           RULES_ALREADY_LOADED = 3};

    try
    {
        convertCMzrPtrToMzrPtr( handle )->attachFileName( std::string(fileName) );
        return SUCCESS;
    }
    catch(utl::dom::noDocumentParsedXcpt xcpt)
   {
       xcpt.warn();
       return DOCUMENT_UNPARSABLE;
   }
   catch(utl::modelAlreadyLoadedXcpt xcpt)
   {
       xcpt.warn();
       return RULES_ALREADY_LOADED;
   }
   catch(utl::xcpt xcpt)
   {
       xcpt.warn();
       return UNKNOWN_ERROR;
   }
    catch(...)
    {
        return UNKNOWN_ERROR;
    }

}

int attachRulesStringToMoleculizerObject( moleculizer* handle, char* rulesCstring)
{
// This function takes a string containing a rules definition of the file and loads 
// it into the provided moleculizer object.

   enum LOCAL_ERROR_TYPE { SUCCESS = 0,
                           UNKNOWN_ERROR=1,
                           DOCUMENT_UNPARSABLE = 2,
                           RULES_ALREADY_LOADED = 3};

   std::string rules( rulesCstring );
   try
   {
       convertCMzrPtrToMzrPtr( handle )->attachString( rules );
       return SUCCESS;
   }
   catch(utl::dom::noDocumentParsedXcpt xcpt)
   {
       xcpt.warn();
       return DOCUMENT_UNPARSABLE;
   }
   catch(utl::modelAlreadyLoadedXcpt xcpt)
   {
       xcpt.warn();
       return RULES_ALREADY_LOADED;
   }
   catch(utl::xcpt xcpt)
   {
       xcpt.warn();
       return UNKNOWN_ERROR;
   }
   catch(...)
   {
       return UNKNOWN_ERROR;
   }

}

int getReactionsBetween(moleculizer* handle, char* cStrSpeciesName1, char* cStrSpeciesName2, reaction*** ptrReactionPtrArray, int* numReactions)
{
    // This function takes two species names, and provides a moleculizer* and a pointer to an array of
    // reaction*'s, and uses moleculizer to determine if there are any reactions between the two species.
    // If so, they are placed in reactionPtrArray and their number is placed in numReactions.

    enum LOCAL_ERROR_TYPE { SUCCESS=0,
                            UNKNOWN_ERROR = 1,
                            ILLEGAL_SPECIES_NAME=2};
    
    std::string speciesName1( cStrSpeciesName1 );
    std::string speciesName2( cStrSpeciesName2 );

    mzr::moleculizer* moleculizerPtr = convertCMzrPtrToMzrPtr( handle );
 
    const mzr::mzrSpecies* species1;
    const mzr::mzrSpecies* species2;

    try
    {
        species1 = moleculizerPtr->getSpeciesWithName( speciesName1 );
        species2 = moleculizerPtr->getSpeciesWithName( speciesName2 );
    }
    catch(mzr::IllegalNameXcpt e)
    {
        return ILLEGAL_SPECIES_NAME; // Species Name is illegal.
    }
    catch(...)
    {
        return UNKNOWN_ERROR;
    }

    std::vector<const mzr::mzrReaction*> reactionsBetweenContainer;
    reactionsBetweenContainer.reserve(10);

    moleculizerPtr->findReactionWithSubstrates( species1, species2, reactionsBetweenContainer);

    if( reactionsBetweenContainer.empty() )
    {
        *ptrReactionPtrArray = NULL;
        *numReactions = 0;
        return SUCCESS;
    }

    typedef reaction* reactionPtr;
    *ptrReactionPtrArray = new reactionPtr[ reactionsBetweenContainer.size() ];

    int size = reactionsBetweenContainer.size();
    *numReactions = size;

    for(unsigned int rxnIndex = 0; rxnIndex != reactionsBetweenContainer.size(); ++rxnIndex)
    {
        (*ptrReactionPtrArray)[ rxnIndex ] = createNewCRxnFromMzrReaction( handle, reactionsBetweenContainer[rxnIndex]);
    }

    return SUCCESS;

}

int getUnaryReactions(moleculizer* handle, char* speciesName, reaction*** ptrReactionPtrArray, int* numReactions)
{
// Same as getBinaryReactions, but with only one species.

    enum LOCAL_ERROR_TYPE { SUCCESS = 0,
                            UNKNOWN_ERROR = 1,
                            ILLEGAL_SPECIES_NAME = 2};

    try
    {
        mzr::moleculizer* moleculizerPtr = convertCMzrPtrToMzrPtr( handle );
        const mzr::mzrSpecies* pSpecies;
    
        pSpecies = moleculizerPtr->getSpeciesWithName( std::string(speciesName) );

        // Just a guess so as to preserve speed when a few reactions are returned;
        std::vector<const mzr::mzrReaction*> unaryReactionContainer;
        unaryReactionContainer.reserve(10);

        moleculizerPtr->findReactionWithSubstrates( pSpecies, unaryReactionContainer);

        if ( unaryReactionContainer.empty() )
        {
            *ptrReactionPtrArray = NULL;
            *numReactions = 0;
            return SUCCESS;
        }

        *numReactions = unaryReactionContainer.size();

        typedef reaction* reactionPtr;        
        *ptrReactionPtrArray = new reactionPtr[ unaryReactionContainer.size() ];

        for(unsigned int rxnIndex = 0; rxnIndex != unaryReactionContainer.size(); ++rxnIndex)
        {
            (*ptrReactionPtrArray)[rxnIndex] = createNewCRxnFromMzrReaction( handle, unaryReactionContainer[ rxnIndex ] );
        }

        return SUCCESS;
    }
    catch( mzr::IllegalNameXcpt )
    {
        return ILLEGAL_SPECIES_NAME;
    }
    catch(utl::xcpt e)
    {
        e.warn();
        return UNKNOWN_ERROR;
    }
    catch(...)
    {
        return UNKNOWN_ERROR;
    }
}

void freeReactionArray( reaction** pRxnArray, unsigned int numElements)
{
    for( unsigned int num = 0; num != numElements; ++num)
    {
        // This will free numElements in the array (all of them);
        freeReaction( pRxnArray[num] );
    }

    // This will delete the allocated memory for the pointers in the array themselves.
    delete [] pRxnArray;
}

void freeSpeciesArray( species** speciesArray, unsigned int numElements)
{
    for(unsigned int num = 0; num != numElements; ++num)
    {
        // Releases the memory each of the pointers points to.
        freeSpecies( speciesArray[num] );
    }

    // Releases the pointers themselves.
    delete [] speciesArray;
}

void freeReaction( reaction* pRxn)
{
    // This frees up all the memory associated with each of the pointers.
    freeSpeciesArray( pRxn->reactantVector, pRxn->numberReactants);
    freeSpeciesArray( pRxn->productVector, pRxn->numberProducts);

    delete pRxn->rate;
}

void freeSpecies( species* pSpecies)
{
    // I am told this is correct under all cases.  
    delete [] pSpecies->name;

    delete pSpecies->mass;
    delete pSpecies->radiusSet;
    delete pSpecies->radius;
    delete pSpecies->diffusionCoeffSet;
    delete pSpecies->diffusionCoeff;

}

void constructCSpeciesArrayFromSpeciesMap(moleculizer* cMzrPtr, 
                                           const mzr::mzrReaction::multMap& speciesMap, 
                                           species*** ptrSpeciesPtrArray, int& numberInList)
{
    int numberSpeciesInMap = getSpeciesMapSize(speciesMap);
    
    typedef species* SpeciesPtr;
    *ptrSpeciesPtrArray = new SpeciesPtr[numberSpeciesInMap] ;

    int counter = 0;
    BOOST_FOREACH( const mzr::mzrReaction::multMap::value_type& value, speciesMap)
    {
        int specCounter = 0;
        while(specCounter++ < value.second)
        {
            // This is the index-th species*
            (*ptrSpeciesPtrArray)[counter++] = createNewCSpeciesFromMzrSpecies(cMzrPtr, value.first);
        }
    }
    
    numberInList = numberSpeciesInMap;

    return;
}


reaction* createNewCRxnFromMzrReaction( moleculizer* cMzrPtr, const mzr::mzrReaction* pMzrReaction)
{
    reaction* theNewRxn = new reaction;
    theNewRxn->rate = new double;

    *(theNewRxn->rate) = pMzrReaction->getRate();

    constructCSpeciesArrayFromSpeciesMap( cMzrPtr, pMzrReaction->getReactants(),&theNewRxn->reactantVector, theNewRxn->numberReactants); 
    constructCSpeciesArrayFromSpeciesMap( cMzrPtr, pMzrReaction->getProducts(), &theNewRxn->productVector, theNewRxn->numberProducts);

    return theNewRxn;

}



int convertUserNameToSpeciesKey( moleculizer* handle, char* theUserName, char* correspondingSpeciesKey)
{

    enum LOCAL_ERROR_TYPE { SUCCESS = 0,
                            NONEXISTANT_USER_NAME = 1,
                            UNKNOWN_ERROR = 2};

    std::string userName( theUserName );

    try
    {
        std::string speciesKey = convertCMzrPtrToMzrPtr(handle)->convertUserNameToGeneratedName( userName );

        // Allocate enough memory for the name as well as a '\0'.  
        correspondingSpeciesKey = new char[ speciesKey.size() + 1];

        // Copy it in.
        strcpy(correspondingSpeciesKey, speciesKey.c_str() );

        assert( strcmp( speciesKey.c_str(), correspondingSpeciesKey) == 0);

        return SUCCESS;

    }
    catch( mzr::unknownUserNameXcpt e)
    {
        e.warn();
        return NONEXISTANT_USER_NAME;
    }
    catch(utl::xcpt e)
    {
        e.warn();
        return UNKNOWN_ERROR;
    }
    catch(...)
    {
        return UNKNOWN_ERROR;
    }
}





int getSpeciesMapSize( const mzr::mzrReaction::multMap& speciesMap)
{
    int size = 0;
    BOOST_FOREACH(const mzr::mzrReaction::multMap::value_type& value, speciesMap)
    {
        size += value.second;
    }
    return size;
}

inline mzr::moleculizer* convertCMzrPtrToMzrPtr(moleculizer* cMzrPtr)
{
    return (mzr::moleculizer*) cMzrPtr->mzrObject;
}


/* The moleculizer here is so that we can look up species SpecificData such as the radius and whatnot. */
species* createNewCSpeciesFromMzrSpecies( moleculizer* cMzrPtr, const mzr::mzrSpecies* pMzrSpecies)
{
    species* newSpecies = new species;

//     newSpecies->name = NULL;
//     newSpecies->mass = NULL;
//     newSpecies->radiusSet = NULL;
//     newSpecies->radius = NULL;
//     newSpecies->diffusionCoeffSet = NULL;
//     newSpecies->diffusionCoeff = NULL;

//     // COPY OVER THE NAME.
//     //
//     std::string speciesKey( pMzrSpecies->getName() );            // Create a string with the name.
//     newSpecies->name = new char[speciesKey.size() + 1];          // Create a new char buffer in newSpecies->name.
//     strcpy( newSpecies->name, speciesKey.c_str() );              // Copy the string into the char buffer.

// #ifdef DEBUG
//     // Double check that the buffers are all working out all right.
//     assert( strcmp( newSpecies->name, speciesKey.c_str() ) == 0 );
// #endif
    
//     // Copy over the mass
//     // 
//     newSpecies->mass = new double;                    // Allocate new mem for a double in mass.
//     *(newSpecies->mass) = pMzrSpecies->getWeight();   // Put the weight there.

//     // Set the radius and the diffusion coefficient
//     newSpecies->radiusSet = new int;                  // Allocate memory for the radisSet.
//     newSpecies->radius = new double;                  // Ditto for the radius.

//     newSpecies->diffusionCoeffSet = new int;          // Allocate memory for the diffusionCoeffSet.
//     newSpecies->diffusionCoeff = new double;          // ... and for the diffusionCoeff.

//     // Attempt to get the radius from the moleculizer handle and set both radiusSet and radius.
//     // If it doesn't work (the radius has not been recorded), record the radiusSet as unset.
//     try
//     {
//         //double theRadiusFromMzr = convertCMzrPtrToMzrPtr(cMzrPtr)->getRadiusForSpecies(pMzrSpecies);
//         //*(newSpecies->radiusSet) = 1;
//         //*(newSpecies->radius) = theRadiusFromMzr;
            
//     }
//     catch(utl::xcpt e)
//     {
//         //*(newSpecies->radiusSet) = 0;
//     }

//     // Attempt to get the kD from moleculizer and use it to set diffusionCoeff.
//     // If it doesn't work, record the diffusionCoeffSet as false.
//     try
//     {
// //         double theDiffusionCoeffFromMzr = convertCMzrPtrToMzrPtr(cMzrPtr)->getKDForSpecies(pMzrSpecies);
// //         *(newSpecies->diffusionCoeffSet) = 1;
// //         *(newSpecies->diffusionCoeff) = theDiffusionCoeffFromMzr;

//     }
//     catch( utl::xcpt e)
//     {
//         *(newSpecies->diffusionCoeffSet) = 0;
//     }
    
    // Return the pointer to the newly created and instantiated newSpecies.
    return newSpecies;
}
