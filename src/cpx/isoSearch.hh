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
// Contact information:
//   Nathan Addy, Research Assistant     Voice: 510-981-8748
//   The Molecular Sciences Institute    Email: addy@molsci.org  
//   2168 Shattuck Ave.                  
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////


#ifndef CPX_ISOSEARCH_HH
#define CPX_ISOSEARCH_HH

#include <iostream>
#include "cpx/plexIso.hh"
#include "cpx/cpxExceptions.hh"
#include "nauty/nauty.h"

namespace cpx
{
    
    template <class plexT>
    class isoSearch
    {

    public:
        isoSearch(const plexT& rLeftPlex,
                  const plexT& rRightPlex) :
            rLeft(rLeftPlex),
            rRight(rRightPlex)
        {
            enforcePlexIsSimpleGraph(rLeft);
            enforcePlexIsSimpleGraph(rRight);
        }


        virtual
        ~isoSearch() 
        {}

        virtual void
        onSuccess(const plexIso& rIso) const
        {}

        bool findIso(void) const;


        // This should be rewritten in a more-better way.  It is all
        // that remains of the old, recursive graph checking code
        // (this can potentially be much, much slower for certain
        // types of complexes (typically highly symmetrical graphs).
        // Nevertheless it is used by the omni code in the ftr unit,

        bool
        findInjection(void) const;

    private:

        typedef typename plexT::molType molType;

        const plexT& rLeft;
        const plexT& rRight;

        bool
        nautyIsomorphismCheck(const plexT& rLeftPlex,
                              const plexT& rRightPlex) const;

        void createGraphFromPlex(const plexT& aPlex, 
                                 graph* theGraph, 
                                 std::vector<int>& labelingMap,
                                 std::vector<int>& ptnMap,
                                 int m) const;

        void
        createMappingBetweenIsomorphicPlexes(plexIso& isoMapping,
                                             const plexT& leftPlex,  std::vector<int>& leftCanonicalLabeling,
                                             const plexT& rightPlex, std::vector<int>& rightCanonicalLabeling) const;

        void enforcePlexIsSimpleGraph(const plexT& plex) const 
            throw(plexIsNotSimpleGraphXcpt);

        bool 
        compareColoringPartitions(const plexT& plex1, const std::vector<int>& labelingMap1, const std::vector<int>& partitionMap1,
                                  const plexT& plex2, const std::vector<int>& labelingMap2, const std::vector<int>& partitionMap2) const;

        void 
        invertMapping( const std::vector<int>& mapToInvert, std::vector<int>& invertedMapping) const;

        // This is used in the findInjection routine.
        // Determines if rCurrentIso can be extended over all the bindings
        // starting at leftBindingIndex in the left plex.  This is the basic
        // recursive step in the process of finding an injection or isomorphism.
        bool
        mapRestBindings(int leftBindingIndex,
                        const plexIso& rCurrentIso) const;

        void
        createPermutationFromRelabeling(const std::vector<int>& labelingMap,
                                        std::vector<int>& permutationMap) const;

        class IndexArraySorter
        {
            const plexT& thePlex;
        public:
            IndexArraySorter( const plexT& plex)
                :
                thePlex(plex)
            {}

            bool operator()(int i, int j)
            {
                return thePlex.mols[i]->getName() < thePlex.mols[j]->getName();
            }
        };

        class bindingIsEqual
        {
            binding firstBinding;
            binding secondBinding;

        public:
            bindingIsEqual( const binding& bnd)
                :
                firstBinding(bnd.first, bnd.second),
                secondBinding( bnd.second, bnd.first)
            {}

            bool operator()(const binding& bnd)
            {
                return (firstBinding == bnd || 
                        secondBinding == bnd);
            }

        private:

        };

    };
}

#include "isoSearchImpl.hh"

#endif
