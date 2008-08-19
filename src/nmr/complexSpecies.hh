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
//   Nathan Addy, Research Assistant     Voice: 510-981-8748
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
#include "namedMolecule.hh"
#include "nmrExceptions.hh"

#include "utl/macros.hh"

#include <iostream>
#include <sstream>
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
        DECLARE_TYPE( std::vector<MinimalMolSharedPtr>, MolList);
        DECLARE_TYPE( int, BndNdx);
        DECLARE_TYPE( MinimalMol::BindingSite,  BindingSite);
        DECLARE_TYPE( MinimalMol::ModificationList, ModificationList)
        typedef std::pair<MolNdx, BndNdx> __HalfBinding;
        DECLARE_TYPE(__HalfBinding, HalfBinding);
        typedef std::pair<HalfBinding, HalfBinding> __Binding;
        DECLARE_TYPE(  __Binding, Binding);
        DECLARE_TYPE( std::vector<Binding>, BindingList)

    public:

        ComplexSpecies();
        ComplexSpecies(ComplexSpeciesCref aComplexSpecies);
        ComplexSpecies(ComplexOutputStateCref aComplexOutputState);

        ~ComplexSpecies()
        {}

        ComplexSpecies& operator=(const ComplexSpecies& crefComplexSpecies);

        void addMolToComplex(MinimalMolSharedPtr ptrMol, AliasCref anAlias) 
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

        std::string repr() const;

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


std::ostream&
operator<<(std::ostream& stream, nmr::ComplexSpeciesCref aComplexSpecies);

#endif

