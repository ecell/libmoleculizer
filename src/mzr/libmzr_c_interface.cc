/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
**
**        This file is part of Libmoleculizer
**
**        Copyright (C) 2001-2009 The Molecular Sciences Institute.
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
#include "mzr/spatialExtrapolationFunctions.hh"
#include <boost/foreach.hpp>

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

void expandNetwork( moleculizer* handle)
{
    mzr::moleculizer* underlyingMoleculizerObject = convertCMzrPtrToMzrPtr( handle );
    underlyingMoleculizerObject->generateCompleteNetwork();
}

int setRateExtrapolation( moleculizer* handle, int extrapolation)
{
    mzr::moleculizer* underlyingMoleculizerObject = convertCMzrPtrToMzrPtr( handle );

    underlyingMoleculizerObject->setRateExtrapolation( extrapolation );
    
    return 0;
}

void freeMoleculizerObject( moleculizer* handle)
{
    // This function frees the c-style moleculizer handle.
    delete convertCMzrPtrToMzrPtr(handle);
    delete handle;
}


int loadRulesFile(moleculizer* handle, char* fileName)
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

int loadRulesString( moleculizer* handle, char* rulesCstring)
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

int getNumberOfSpecies(moleculizer* handle)
{
    mzr::moleculizer* moleculizerPtr = convertCMzrPtrToMzrPtr( handle );
    int numberOfSpecies = moleculizerPtr->getTotalNumberSpecies();

    return numberOfSpecies;
}

int getNumberOfReactions(moleculizer* handle)
{
    mzr::moleculizer* moleculizerPtr = convertCMzrPtrToMzrPtr( handle );
    int numberOfReactions = moleculizerPtr->getTotalNumberReactions();

    return numberOfReactions;
}

int getAllStreamSpecies(moleculizer* handle, char* cStrStreamName, species*** pSpeciesArray, int* numberSpecies)
{
    enum LOCAL_ERROR_TYPE { SUCCESS = 0,
                            UNKNOWN_ERROR = 1,
                            STREAM_DOES_NOT_EXIST = 3 };

    mzr::moleculizer* moleculizerPtr = convertCMzrPtrToMzrPtr( handle );

    std::string streamName( cStrStreamName );
    std::vector<const mzr::mzrSpecies*> speciesVector;

    moleculizerPtr->getSpeciesInSpeciesStream( streamName, speciesVector );

    *numberSpecies = moleculizerPtr->getTotalNumberSpecies();
    *pSpeciesArray = new species* [ moleculizerPtr->getTotalNumberSpecies() ];


    int index = 0;
    BOOST_FOREACH( const mzr::mzrSpecies* pSpecies, speciesVector)
    {
        (*pSpeciesArray)[index] = createNewCSpeciesFromMzrSpecies( handle, pSpecies);
    }

    return SUCCESS;
}

int getAllSpecies(moleculizer* handle, species*** pSpeciesArray, int* numberSpecies)
{
    enum LOCAL_ERROR_TYPE { SUCCESS = 0,
                            UNKNOWN_ERROR = 1,
                            STREAM_DOES_NOT_EXIST = 3 };

    mzr::moleculizer* moleculizerPtr = convertCMzrPtrToMzrPtr( handle );

    *numberSpecies = moleculizerPtr->getTotalNumberSpecies();
    *pSpeciesArray = new species* [ *numberSpecies ];

    int index = 0;
    BOOST_FOREACH( const mzr::moleculizer::SpeciesCatalog::value_type& refVT, moleculizerPtr->getSpeciesCatalog())
    {
        (*pSpeciesArray)[index++] = createNewCSpeciesFromMzrSpecies( handle, refVT.second );
    }

    return SUCCESS;
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
    reactionsBetweenContainer.reserve(4);
    
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
    delete pSpecies->radius;
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
    
    newSpecies->name = NULL;
    newSpecies->mass = NULL;
    newSpecies->radius = NULL;
    newSpecies->diffusionCoeff = NULL;
    
    // COPY OVER THE NAME.
    std::string speciesKey( pMzrSpecies->getName() );            // Create a string with the name.
    newSpecies->name = new char[speciesKey.size() + 1];          // Create a new char buffer in newSpecies->name.
    strcpy( newSpecies->name, speciesKey.c_str() );              // Copy the string into the char buffer.
    
    // Copy over the mass
    // 
    newSpecies->mass = new double;                    // Allocate new mem for a double in mass.
    *(newSpecies->mass) = pMzrSpecies->getWeight();   // Put the weight there.
    
    // Set the radius and the diffusion coefficient
    newSpecies->radius = new double;      
    newSpecies->diffusionCoeff = new double;          // ... and for the diffusionCoeff.
    
    // Extrapolate the radius and diffusion coefficient from moleculizer.
    *newSpecies->radius = mzr::extrapolateMolecularRadius( pMzrSpecies );
    *newSpecies->diffusionCoeff = mzr::getDiffusionCoeffFromSpecies( pMzrSpecies );
    
    // Return the pointer to the newly created and instantiated newSpecies.
    return newSpecies;
}
