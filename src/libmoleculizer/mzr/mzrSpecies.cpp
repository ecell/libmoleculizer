//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2009 The Molecular Sciences Institute.
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

#include "mzr/mzrSpecies.hpp"
#include "mzr/mzrUnit.hpp"
#include "mzr/unitsMgr.hpp"
#include "mzr/moleculizer.hpp"
#include "mzr/mzrSpeciesDumpable.hpp"
#include "plex/mzrPlexSpecies.hpp"
#include <iostream>

namespace mzr
{
    plx::mzrPlexSpecies*
    mzrSpecies::getComplexSpecies()
    {
	plx::mzrPlexSpecies* ptrComplexSpecies = dynamic_cast<plx::mzrPlexSpecies*>(this);
	return ptrComplexSpecies;
    }

    void
    mzrSpecies::expandReactionNetwork()
    {
        this->expandReactionNetwork( 0 );
    }
    
    void
    mzrSpecies::expandReactionNetwork( unsigned int depth )
    {
        ensureNotified( depth );
    }
    
    void
    mzrSpecies::setGenerateDepth( unsigned int i )
    {
        mzrSpecies::generateDepth = i;
    }
    
    
    unsigned int
    mzrSpecies::getGenerateDepth()
    {
        return mzrSpecies::generateDepth;
    }

    unsigned int mzrSpecies::generateDepth = 1;
}
