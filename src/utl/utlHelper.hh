/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2008 Nathan Addy
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
/////////////////////////////////////////////////////////////////////////////

#ifndef UTLHELPER_HH
#define UTLHELPER_HH

#include <utility>

namespace utl
{
    namespace aux
    {

        // I want to create two functional objects: 
        // 1. One which is given a pair type and an int 
        //    (getNthPositionFromPair< std::pair<string*, mzrSpecies*>, 2>
        //    and returns the object in that position

        // 2.  Something that takes a * and returns the dereferenced value.

//         template <typename Tunsigned int>
//         class getNthPositionFromPair
//         {
//         public:
//             template
//         }

        class getFirstFromPair
        {

            template <typename T, typename V>
            const T& operator()(const std::pair<T, V>& aPair)
            {
                
            }
        };

        class ptrDereferenceSorter
        {
        public:
            template <typename T>
            bool operator()(const T* a, const T* b) const
            {
                return (*a) < (*b);
            }
        };
        
    }
}

#endif
