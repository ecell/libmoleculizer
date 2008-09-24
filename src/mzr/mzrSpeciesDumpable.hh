//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                                                                          
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

#ifndef MZR_MZRSPECIESDUMPABLE_H
#define MZR_MZRSPECIESDUMPABLE_H

#include "utl/dom.hh"
#include "fnd/varDumpable.hh"
#include "fnd/querySpeciesDumpable.hh"
#include "mzr/mzrSpecies.hh"
#include "mzr/mzrStream.hh"

namespace mzr
{
// These templates add routines to emit output for dumpables
// in a simulation state dump.  This keeps xmlpp stuff out of the
// main templates.
template<class mzrSpeciesType>
class singleSpeciesDumpable :
public fnd::varDumpable<mzrSpeciesType,
fnd::basicDumpable::dumpArg>,
public mzrSpeciesStream
{
public:
singleSpeciesDumpable(const std::string& rName,
const mzrSpeciesType* pSpeciesToDump) :
fnd::varDumpable<mzrSpeciesType,
fnd::basicDumpable::dumpArg>(rName,
pSpeciesToDump)
{}

~singleSpeciesDumpable(void)
{}

void
insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
throw(std::exception);

void
insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
throw(std::exception);
};

template<class mzrSpeciesType>
class multiSpeciesDumpable :
public fnd::multiSpeciesDumpable<mzrSpeciesType,
fnd::basicDumpable::dumpArg>,
public mzrSpeciesStream
{
public:
multiSpeciesDumpable(const std::string& rName) :
fnd::multiSpeciesDumpable<mzrSpeciesType,
fnd::basicDumpable::dumpArg>(rName)
{}

~multiSpeciesDumpable(void)
{}

void
insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
throw(std::exception);

void
insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
throw(std::exception);
};

template<class mzrSpeciesType>
class querySpeciesDumpable :
public fnd::querySpeciesDumpable<mzrSpeciesType,
fnd::basicDumpable::dumpArg>,
public mzrSpeciesStream
{
public:
querySpeciesDumpable(const std::string& rName,
fnd::query<mzrSpeciesType>& rQuery) :
fnd::querySpeciesDumpable<mzrSpeciesType,
fnd::basicDumpable::dumpArg>(rName,
rQuery)
{}

~querySpeciesDumpable(void)
{}

void
insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
throw(std::exception);

void
insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
throw(std::exception);
};
}

#include "mzr/mzrSpeciesDumpableImpl.hh"

#endif // MZR_MZRSPECIESDUMPABLE_H
