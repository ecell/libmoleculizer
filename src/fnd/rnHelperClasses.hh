//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2008 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Nathan Addy, Scientific Programmer, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//


#ifndef RNHELPERCLASSES
#define RNHELPERCLASSES

#include <iostream>

namespace fnd
{
namespace aux
{
template <typename T>
class compareByPtrValue
{
public:
bool operator()(const T* const a, const T* const b) const
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
