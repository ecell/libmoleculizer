#include <vector>
#include "modMol.hh"
#include "mzr/debug.hh"
#include "utl/string.hh"
#include "binding.hh"

// Can't believe I have to redo all this... 
//       --NJA

namespace cpx
{
    template <class plexFamilyT>
    void 
    plexSpeciesMixin<plexFamilyT>::
    createComplexRepresentation( ComplexRepresentation& aComplexRepresentation ) const
    {
        // Is there even a way to access this???
        // static std::vector<std::string> allLegalModifications = 

        // There is probably some way to clear() a ComplexRepresentation, and we should 
        // do that here.  
        // aComplexRepresentation.clear();

        const std::vector<typename plexFamilyT::molType*>& rMols = rFamily.getParadigm().mols;
        const std::vector<cpx::binding>& rBindings  = rFamily.getParadigm().bindings;

        // Add the mols to the complex representation.
        for(unsigned int molNdx = 0;
            molNdx != rMols.size();
            ++molNdx)
        {
            typename plexFamilyT::molType* pMol = rMols[molNdx];
            
            complexspecies::SimpleMol aMol( pMol->getName() );

            // Create all the binding sites to this mol.
            for( typename plexFamilyT::molType::bindingSiteIterator iter = pMol->getBindingSitesBegin();
                 iter != pMol->getBindingSitesEnd();
                 ++iter)
            {
                aMol.addNewBindingSite( (*iter).getName() );
            }

            // This line seems to confuse the bejeezus out of emacs, and I don't know why.
            const cpx::modMol<typename plexFamilyT::molType>* aModMol = dynamic_cast<const cpx::modMol<typename plexFamilyT::molType>* >(pMol);

            if(aModMol)
            {
                // Get the externalized state....
                const cpx::modMolState& nuMolParam = aModMol->externState( molParams[molNdx] );

                 if( nuMolParam.size() !=aModMol->modSiteNames.size() )
                 {
                     // Throw something better...
                     throw 666;
                 }
                


                for(unsigned int ndx = 0;
                    ndx != aModMol->modSiteNames.size();
                    ++ndx)
                {
                    
                    aMol.addNewModificationSite( aModMol->modSiteNames[ndx],
                                                 nuMolParam[ndx]->getName() );
                }

            }

            aComplexRepresentation.addMolToComplex( aMol, utl::stringify(molNdx) );            

        }
    

     // Add the bindings to the complex representation.
     for(std::vector<cpx::binding>::const_iterator iter = rBindings.begin();
         iter != rBindings.end();
         ++iter)
     {
         int mol1_Ndx = (*iter).leftSite().molNdx();
         int mol2_Ndx = (*iter).rightSite().molNdx();

         int mol1_siteNdx = (*iter).leftSite().siteNdx();        
         int mol2_siteNdx = (*iter).rightSite().siteNdx();


         std::string site1Name = (*rMols[mol1_Ndx])[mol1_siteNdx].getName();
         std::string site2Name = (*rMols[mol2_Ndx])[mol2_siteNdx].getName();


         aComplexRepresentation.addBindingToComplex( utl::stringify(mol1_Ndx),
                                                     site1Name,
                                                     utl::stringify(mol2_Ndx),
                                                     site2Name );
     }   
     
        return;
    }


    template <class plexFamilyT>
    std::string
    plexSpeciesMixin<plexFamilyT>::getCanonicalName(void) const
    {
        // static complexspecies::basicNameAssembler<complexspecies::SimpleMol> defaultNameAssembler;
        static complexspecies::MangledNameAssembler<complexspecies::SimpleMol> defaultNameAssembler;

        // complexspecies::readableNameAssembler<complexspecies::SimpleMol> defaultNameAssembler;
        return getCanonicalName( &defaultNameAssembler );
    }

    template <class plexFamilyT>
    std::string
    plexSpeciesMixin<plexFamilyT>::
    getCanonicalName( const complexspecies::NameAssembler<complexspecies::SimpleMol>* const ptrNameAssembler) const
    {
        complexspecies::ComplexSpecies<complexspecies::SimpleMol> aComplexSpecies;
        createComplexRepresentation(aComplexSpecies);
        string theName( ptrNameAssembler->createCanonicalName(aComplexSpecies) );

        return theName;
    }

}
