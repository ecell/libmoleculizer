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
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001
//
// Modifing Authors:
//
//

#include "mzr/mzrEltName.hh"

// Aren't these all the same? Could they go into mzrSpeciesStream?

namespace mzr
{
template<class mzrSpeciesType>
void
singleSpeciesDumpable<mzrSpeciesType>::
insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
throw(std::exception)
{
xmlpp::Element* pTaggedSpeciesStreamRefElt
= pParent->add_child(eltName::taggedSpeciesStreamRef);

pTaggedSpeciesStreamRefElt
->set_attribute(eltName::taggedSpeciesStreamRef_nameAttr,
this->getName());
}

template<class mzrSpeciesType>
void
singleSpeciesDumpable<mzrSpeciesType>::
insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
throw(std::exception)
{
xmlpp::Element* pSpeciesRefElt
= pParentElt->add_child(eltName::taggedSpeciesRef);

pSpeciesRefElt->set_attribute(eltName::taggedSpeciesRef_tagAttr,
this->getVar()->getTag());
}

template<class mzrSpeciesType>
void
multiSpeciesDumpable<mzrSpeciesType>::
insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
throw(std::exception)
{
xmlpp::Element* pTaggedSpeciesStreamRefElt
= pParent->add_child(eltName::taggedSpeciesStreamRef);

pTaggedSpeciesStreamRefElt
->set_attribute(eltName::taggedSpeciesStreamRef_nameAttr,
this->getName());
}

template<class mzrSpeciesType>
class insertMultiSpeciesDumpableTag :
public std::unary_function<mzrSpeciesType*, void>
{
xmlpp::Element* pParentElt;

public:
insertMultiSpeciesDumpableTag(xmlpp::Element* pParentElement) :
pParentElt(pParentElement)
{}

void
operator()(const mzrSpeciesType* pSpecies) const
throw(std::exception)
{
xmlpp::Element* pTaggedSpeciesRefElt
= pParentElt->add_child(eltName::taggedSpeciesRef);

pTaggedSpeciesRefElt->set_attribute(eltName::taggedSpeciesRef_tagAttr,
pSpecies->getTag());
}
};

template<class mzrSpeciesType>
void
multiSpeciesDumpable<mzrSpeciesType>::
insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
throw(std::exception)
{
std::for_each(this->dumpedSpecies.begin(),
this->dumpedSpecies.end(),
insertMultiSpeciesDumpableTag<mzrSpeciesType>(pParentElt));
}

template<class mzrSpeciesType>
void
querySpeciesDumpable<mzrSpeciesType>::
insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
throw(std::exception)
{
xmlpp::Element* pTaggedSpeciesStreamRefElt
= pParent->add_child(eltName::taggedSpeciesStreamRef);

pTaggedSpeciesStreamRefElt
->set_attribute(eltName::taggedSpeciesStreamRef_nameAttr,
this->getName());
}

template<class mzrSpeciesType>
class insertQuerySpeciesDumpableTag :
public std::unary_function<mzrSpeciesType*, void>
{
xmlpp::Element* pParentElt;

public:
insertQuerySpeciesDumpableTag(xmlpp::Element* pParentElement) :
pParentElt(pParentElement)
{}

void
operator()(const mzrSpeciesType* pSpecies) const
throw(std::exception)
{
xmlpp::Element* pTaggedSpeciesRefElt
= pParentElt->add_child(eltName::taggedSpeciesRef);

pTaggedSpeciesRefElt->set_attribute(eltName::taggedSpeciesRef_tagAttr,
pSpecies->getTag());
}
};

template<class mzrSpeciesType>
void
querySpeciesDumpable<mzrSpeciesType>::
insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
throw(std::exception)
{
std::for_each(this->dumpedSpecies.begin(),
this->dumpedSpecies.end(),
insertQuerySpeciesDumpableTag<mzrSpeciesType>(pParentElt));
}
}
