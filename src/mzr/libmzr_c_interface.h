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


#ifdef __cplusplus
extern "C" {
#endif
    
    typedef struct moleculizer_handle
    {
        void* mzrObject;
        
    } moleculizer;
    
    typedef struct species_type
    {
        char* name;
        double* mass;  /* These are in daltons */
        
        int* radiusSet;  
        double* radius; 
        
        int* diffusionCoeffSet;
        double* diffusionCoeff;
        
    } species;
    
    typedef struct reaction_type
    {
        int numberReactants;
        species** reactantVector;
        
        int numberProducts;
        species** productVector;
        
        double* rate;
        
    } reaction;
    
    
    /* These fuctions create a new moleculizer object and release it after being used.
       A call to createNewMoleculizerObject should always be paired with a call to 
       freeMoleculizerObject */
    
    moleculizer* createNewMoleculizerObject();
    
    int setupSpatialExtrapolation(moleculizer* handle);
    int setupNonSpatialExtrapolation(moleculizer* handle, int useMassBasedRateExtrapolation);
    
    void freeMoleculizerObject( moleculizer* handle);
    
    
    
    /* This function should be called with a file that contains an xml file describing the rules of the 
       system. See documentation for a description of the rules. */
    
    int attachRulesFileToMoleculizerObject(moleculizer* handle, char* fileName);
    int attachRulesStringToMoleculizerObject( moleculizer* handle, char* file);
    
    
    /* These two functions can be used with a Species Key (it's canonical string representation) to 
       determine what reactions, if any they participate in. These functions return pointers to reaction 
       arrays, which contain pointers to reactions and species.  Ownership is transfered to the user, 
       who must free them manually. */
    
    int getReactionsBetween(moleculizer* handle, char* speciesName1, char* speciesName2, reaction*** ptrReactionPtrArray, int* numReactions);
    int getUnaryReactions(moleculizer* handle, char* speciesName, reaction*** ptrReactionPtrArray, int* numReactions);
    
    
    
    /* These functions free species and reaction arrays, such as those provided by the getReactionsBetween
       and getUnaryReactions functions.  Usually, a call to freeReaction( reactionArray, *numReactions)
       is sufficient to properly release all this information. */
    
    void freeReactionArray( reaction** pRxnArray, unsigned int numArrayElements);
    void freeSpeciesArray( species** pSpeciesArray, unsigned int numArrayElements);
    
    void freeReaction( reaction* pRxn );
    void freeSpecies( species* pSpecies );
    
    
    /* This function converts a user name to a species key, ie like "X-singleton" to "___1X______" or 
       whatever.
       
       I may have to change this, as I am not certain of the best c-interface here. For the moment, 
       theUserName must be null terminated.  It returns non-zero if in error (usually it means
       the name has not been found.. */
    
    int convertUserNameToSpeciesKey(char* theUserName, char* correspondingSpeciesKey);
    
    
#ifdef __cplusplus
}
#endif
