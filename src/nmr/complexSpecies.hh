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

#ifndef __COMPLEXSPECIES_HH
#define __COMPLEXSPECIES_HH

#include "complexOutputState.hh"
#include "partialTokenList.hh"
#include "permutation.hh"
#include "abstractMol.hh"
#include "nmrExceptions.hh"

#include "utl/macros.hh"

#include <string>
#include <vector>
#include <utility>
#include <iterator>
#include <map>

namespace nmr
{


    // April 20 - I am proceeding to change the Mols that ComplexSpecies have to 
    // boost::shared_ptr<Mol> in the hopes of more efficient memory management 
    // shoon/eventually.
    DECLARE_CLASS( ComplexSpecies );
    class ComplexSpecies
    {

    public:
        DECLARE_TYPE( std::string, Alias);
        DECLARE_TYPE( int, MolNdx );
        DECLARE_TYPE( std::vector<MolSharedPtr>, MolList);
        DECLARE_TYPE( int, BndNdx);
        DECLARE_TYPE( Mol::BindingSite,  BindingSite);
        DECLARE_TYPE( Mol::ModificationList, ModificationList)
        typedef std::pair<std::pair<MolNdx, BndNdx>, std::pair<MolNdx, BndNdx> > __Binding;
        DECLARE_TYPE(  __Binding, Binding);
        DECLARE_TYPE( std::vector<Binding>, BindingList)

    public:

        ComplexSpecies();
        ComplexSpecies(ComplexSpeciesCref aComplexSpecies);
        ComplexSpecies(ComplexOutputStateCref aComplexOutputState);

        ~ComplexSpecies()
        {}

        void addMolToComplex(MolSharedPtr ptrMol, AliasCref anAlias) 
            throw(DuplicateMolAliasXcpt);

        void addBindingToComplex(AliasCref firstMolAlias, 
                                 BindingSiteCref firstMolSiteAlias, 
                                 AliasCref secondMolAlias, 
                                 BindingSiteCref secondMolSiteAlias) throw( MissingMolAliasXcpt, 
                                                                            MissingBindingSiteXcpt );

        unsigned int 
        getNumberOfMolsInComplex() const;

        unsigned int 
        getNumberOfBindingsInComplex() const;
  
        MolListCref 
        getMolList() const;
        
        BindingListCref 
        getBindingList() const;

        MolListRef 
        getMolList();
        
        BindingListRef 
        getBindingList();
        
        void 
        applyPermutationToComplex(const Permutation& aPermutation);

        // This is a partialNameSentence with everything "stringified".  This plexOutputState is 
        // the object passed to a particular nameAssembler.
        void constructOutputState(ComplexOutputState& rOutputState) const;

        void 
        DEBUG_print() const
        {
            std::cout << "Printing theMols" << std::endl;
            for(MolList::const_iterator iter = theMols.begin();
                iter != theMols.end();
                ++iter)
            {
                std::cout << (*iter)->getMolType() << std::endl;
            }
            std::cout << "##############################" << std::endl;

            std::cout << "Printing the Alias to Ndx map" << std::endl;
            for( MolMap::const_iterator iter = theMolAliasToNdxMap.begin();
                 iter != theMolAliasToNdxMap.end();
                 ++iter)
            {
                std::cout << iter->first << '\t' << iter->second << std::endl;
            }
        }

    protected:
        typedef std::map<Alias, MolNdx> _molMap;
        DECLARE_TYPE( _molMap, MolMap);

        typedef MolMap::iterator MolMapIter;
        DECLARE_TYPE(ComplexOutputState::MolTokenStr, MolTokenStr);
        DECLARE_TYPE(ComplexOutputState::BindingTokenStr, BindingTokenStr);
        DECLARE_TYPE(ComplexOutputState::ModificationTokenStr, ModificationTokenStr);

        void constructPartialTokenList(PartialTokenListRef rComplexPartialTokenList) const;
        void sortBinding(BindingRef aBinding);

        MolMap theMolAliasToNdxMap;  
        MolList theMols;
        BindingList theBindings;

    };

}

#endif

