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

#ifndef MZR_SIMTIMEDUMPABLE_H
#define MZR_SIMTIMEDUMPABLE_H

#include "fnd/dumpable.hh"

namespace mzr
{
class moleculizer;

/*! \ingroup dumpGroup
\brief Dumps the simulation time.

Automatically used as the first field by most dump commands. */
class simTimeDumpable :
            public fnd::dumpable<fnd::basicDumpable::dumpArg>
{
    mzr::moleculizer& rMolzer;

public:
    simTimeDumpable( mzr::moleculizer& rMoleculizer );

    virtual void
    doDump( const fnd::basicDumpable::dumpArg& rDumpArg ) const;
};
}

#endif // MZR_SIMTIMEDUMPABLE_H
