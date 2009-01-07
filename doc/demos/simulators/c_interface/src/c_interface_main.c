#include <stdio.h>
#include "mzr/libmzr_c_interface.h"

int main(int argc, char* argv[] )
{
  moleculizer* pMoleculizer = createNewMoleculizerObject();

  if (argc < 3 || strcmp(argv[1], "-f") != 0 )
    {
      printf("No model file was loaded.\n");
      return 1;
    }
  else
    {
      int code = loadRulesFile( pMoleculizer, argv[ 2 ] );
      switch(code)
	{
	case 0:
	  printf("Model loaded successfully.\n");
	  break;
	case 1:
	  printf("Unknown error on load.  Aborting\n");
	  return 1;
	case 2:
	  printf("Document unparsable.  Aborting.\n");
	  return 1;
	case 3:
	  printf("Moleculizer has already loaded rules. Ignoring and continuing.\n");
	  return 1;
	}
    }

    setRateExtrapolation( pMoleculizer, 1 );

    expandNetwork( pMoleculizer );

    printf("Expanded entire network.\n");

    printf("There are %d species and %d reactions.\n", 
           getNumberOfSpecies(pMoleculizer),
           getNumberOfReactions(pMoleculizer) );
           

    freeMoleculizerObject( pMoleculizer );
    return 0;
}


