/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008 The Molecular Sciences Institute
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
/////////////////////////////////////////////////////////////////////////////

#include "utl/stdIncludes.hh"
#include "utl/utility.hh"

#include "modMol.hh"
#include "binding.hh"

// Veni, Vidi, Vici. -- Julius Ceasar 

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
            
            nmr::MinimalMolSharedPtr aMol( new nmr::MinimalMol( pMol->getName() ) );

            // Create all the binding sites to this mol.
            for( typename plexFamilyT::molType::bindingSiteIterator iter = pMol->getBindingSitesBegin();
                 iter != pMol->getBindingSitesEnd();
                 ++iter)
            {
                aMol->addNewBindingSite( (*iter).getName() );
            }

            const cpx::modMol<typename plexFamilyT::molType>* aModMol = 
                dynamic_cast<const cpx::modMol<typename plexFamilyT::molType>* >(pMol);

            if(aModMol)
            {
                // Get the externalized state....
                const cpx::modMolState& nuMolParam = aModMol->externState( molParams[molNdx] );

                 if( nuMolParam.size() !=aModMol->modSiteNames.size() )
                 {
                     throw utl::xcpt("Unknown Error in plexSpeciesMixin::createComplexRepresentation. (key: axjfdek)");
                 }
                
                for(unsigned int ndx = 0;
                    ndx != aModMol->modSiteNames.size();
                    ++ndx)
                {
                    
                    aMol->addNewModificationSite( aModMol->modSiteNames[ndx],
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
        // static nmr::basicNameAssembler defaultNameAssembler;
        static nmr::MangledNameAssembler defaultNameAssembler;
        // nmr::readableNameAssembler defaultNameAssembler;

        return getCanonicalName( &defaultNameAssembler );
    }

    template <class plexFamilyT>
    std::string
    plexSpeciesMixin<plexFamilyT>::
    getCanonicalName( const nmr::NameAssembler* const ptrNameAssembler) const
    {
        nmr::ComplexSpecies aComplexSpecies;
        createComplexRepresentation(aComplexSpecies);
        string theName( ptrNameAssembler->createCanonicalName(aComplexSpecies) );
        
        return theName;
    }

}
