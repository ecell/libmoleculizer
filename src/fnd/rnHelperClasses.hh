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


#ifndef RNHELPERCLASSES
#define RNHELPERCLASSES

namespace fnd
{
    namespace aux
    {

        template <typename T>
        class compareByPtrValue
        {
        public:
            bool operator()(const T* a, const T* b)
            {
                return *a < *b;
            }
        };

        template <typename ListCatalogT>
        class doDeleteStringPtrs
            : public std::unary_function<typename ListCatalogT::value_type, void>
        {
        public:
            void operator()(const typename doDeleteStringPtrs::argument_type& refPairWithString)
            {
                delete refPairWithString.first;
            }
        };

        template <typename ListCatalogT>
        class findEntryWithName
            : public std::unary_function<typename ListCatalogT::value_type, bool>
        {
        public:
            findEntryWithName(const std::string& nameToFind)
                :
                nameMatchTarget( nameToFind )
            {}

            bool operator()(const typename findEntryWithName::argument_type& potentialMatch)
            {
                return (*potentialMatch.first == nameMatchTarget);
            }

        private:
            const std::string& nameMatchTarget;
        };
    }
}


#endif // RNHELPERCLASSES
