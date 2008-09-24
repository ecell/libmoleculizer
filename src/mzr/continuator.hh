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

#ifndef CONTUNUATOR_H
#define CONTUNUATOR_H

#include "utl/dom.hh"
#include "mzr/moleculizer.hh"

namespace mzr
{
class continuator :
public moleculizer
{
public:
// The point of this program is to continue a moleculizer simulation
// given the original moleculizer-input and a moleculizer-state dump.
//
// It parses moleculizer-input as usual, then parses tagged species from
// moleculizer-state.  The populations of the tagged species are used for
// all species, instead of the populations given for explicit-species in
// explicit-species.  The initial simulation time is set to the simulation
// time at which state was dumped.  Thus the simulation picks up where it
// left off.
continuator(int argc,
char** argv,
xmlpp::Document* pMoleculizerInput,
xmlpp::Document* pMoleculizerState)
throw(std::exception);

~continuator(void)
{}
};
}

#endif




