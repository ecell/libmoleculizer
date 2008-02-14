/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001  Walter Lawrence (Larry) Lok.
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
//   Larry Lok, Research Fellow          Voice: 510-981-8740
//   The Molecular Sciences Institute      Fax: 510-647-0699
//   2168 Shattuck Ave.                  Email: lok@molsci.org
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#include "fnd/sensitivityList.hh"
#include "mzr/mzrSpecies.hh"
#include "mzr/createEvent.hh"
#include "mzr/respondReaction.hh"
#include "mzr/mzrUnit.hh"

namespace mzr
{
  // This kind of event is for initializing the simulation;
  // it just instantly creates the given number of the given
  // species.
  createEvent::createEvent(mzrSpecies* pSpecies,
			   int count,
			   mzrUnit& refMzrUnit) :
    pSpeciesToCreate(pSpecies),
    howMany(count),
    rMzrUnit(refMzrUnit)
  {}

  fnd::eventResult
  createEvent::happen(moleculizer& rMolzer)
    throw(std::exception)
  {
    int generateDepth
      = rMzrUnit.getGenerateDepth();
    
    fnd::sensitivityList<mzrReaction> affectedReactions;
    pSpeciesToCreate->update(howMany,
			     affectedReactions,
			     generateDepth);

    for_each(affectedReactions.begin(),
	     affectedReactions.end(),
	     respondReaction(rMolzer));

    return fnd::go;
  }
}
