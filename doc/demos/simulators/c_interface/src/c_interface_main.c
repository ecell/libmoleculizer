#include <stdio.h>
#include "mzr/libmzr_c_interface.h"

int processArgs( int argc, char* argv[], moleculizer* pMzr);

int main(int argc, char* argv[] )
{

    moleculizer* pMoleculizer = createNewMoleculizerObject();

    if (!processArgs( argc, argv, pMoleculizer))
    {
        printf("No model file was loaded...\n");
        return 1;
    }
    setRateExtrapolation( pMoleculizer, 1 );

    expandNetwork( pMoleculizer );

    printf("Expanded entire network.\n");

    printf("There are %d species and %s reactions.\n", 
           getNumberOfSpecies(pMoleculizer),
           getNumberOfReactions(pMoleculizer) );
           

    freeMoleculizerObject( pMoleculizer );
    return 0;
}




int processArgs( int argc, char* argv[], moleculizer* handle)
{
    int loaded = 0;
    
    int i;
    for (i = 1; i < argc; i++)  /* Skip argv[0] (program name). */
    {
        /*
         * Use the 'strcmp' function to compare the argv values
         * to a string of your choice (here, it's the optional
         * argument "-q").  When strcmp returns 0, it means that the
         * two strings are identical.
         */

        if (strcmp(argv[i], "-f") == 0)  /* Process optional arguments. */
        {
            loadRulesFile(handle, argv[i+1]);
            loaded = 1;
        }
    }

    return loaded;
}

