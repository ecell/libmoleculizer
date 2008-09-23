/////////////////////////////////////////////////////////////////////////////
// libComplexSpecies - a library for canonically naming species of protein 
//                     complexes.
// Copyright (C) 2008 The Molecular Sciences Institute
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

#ifndef COMPLEXSPECIESOUTPUTMINIMIZER_HPP
#define COMPLEXSPECIESOUTPUTMINIMIZER_HPP


#include "complexSpecies.hh"
#include "permutation.hh"
#include <vector>
#include <set>


namespace nmr
{
    class ComplexSpeciesOutputMinimizer
    {
    public:

        typedef std::vector<int> ColoringPartition;
        typedef std::vector< std::set<int>* > GraphEdgeList;

        ComplexSpeciesOutputMinimizer()
        {}

        ComplexOutputState
        getMinimalOutputState(ComplexSpeciesCref aComplexSpecies);

        class NonSimpleGraphXcpt : public GeneralNmrXcpt
        {
        public:
            static 
            std::string
            mkMsg(ComplexSpeciesCref aComplexSpecies)
            {
                std::ostringstream oss;
                oss << "Error: the complex species '"
                    << aComplexSpecies 
                    << "' does not represent a simple graph.  Please contact "
                    << "the developer <addy@molsci.org>.";
                return oss.str();
            }
            
            NonSimpleGraphXcpt(ComplexSpeciesCref aComplexSpecies)
                :
                GeneralNmrXcpt(mkMsg(aComplexSpecies))
            {}
        };

    private:
        Permutation 
        calculateMolSortingPermutationForComplex( ComplexSpeciesCref aComplexSpecies);
        
        Permutation 
        calculateCanonicalPermutationForColoredGraph( const GraphEdgeList& refGraphEdgeList,
                                                      const ColoringPartition& refColPart);

        bool checkComplexSpeciesIsSimpleGraph( ComplexSpeciesCref aComplexSpecies);
        void setupDataStructuresForCalculation( ComplexSpeciesRef aComplexSpecies);
        void setupComplexEdgeMap( ComplexSpeciesCref aComplexSpecies);
        void setupComplexColorPartition( ComplexSpeciesCref aComplexSpecies);

        struct MolIndexLessThanCmp : public std::binary_function<int, int, bool>
        {
            DECLARE_TYPE(ComplexSpecies::MolList, MolList);

            MolIndexLessThanCmp(ComplexSpeciesCref aComplexSpeciesForCmp);
            bool operator()(int ndx1, int ndx2);

        private:
            MolListCref theComparisonMolList;
        };  

    private:
        ColoringPartition partitionSpecification;
        GraphEdgeList complexGraphEdgeMap;
        // std::vector<std::set<int>* > complexGraphEdgeMap;
    };
}






#endif
