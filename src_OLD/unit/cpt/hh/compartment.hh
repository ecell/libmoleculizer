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

#ifndef COMPARTMENT_H
#define COMPARTMENT_H

#include <string>
#include "utl/dom.hh"
#include "fnd/physConst.hh"
#include "cpt/cptReaction.hh"

namespace cpt
{
  class compartment
  {
    // The name given in the input file.
    std::string name;

    // The volume given in the input file.
    fnd::stateVar<double,
		  cptReaction> volume;

    // (With reservations) the index of this compartment in its
    // compartmentGraph.
    int index;

  public:
    // For construction by parser.  Does not set the compartment index,
    // which is supposed to be set when the compartment is added to a
    // compartment graph.
    compartment(const std::string& rName,
		double initialVolume) :
      name(rName),
      volume(initialVolume),
      index(-1)
    {}

    const std::string&
    getName(void) const
    {
      return name;
    }

    void
    setIndex(int newIndex)
    {
      index = newIndex;
    }

    int
    getIndex(void) const
    {
      return index;
    }

    void
    updateVolume(double newVolume,
		 fnd::sensitivityList<cptReaction>& rAffectedReactions)
    {
      volume.updateValue(newVolume,
			 rAffectedReactions);
    }

    // This may not be a good strategy.
    void
    sensitizeReactionToVolume(cptReaction& rReaction)
    {
      volume.addSensitive(&rReaction);
    }

    double
    getVolume(void) const
    {
      return volume.getValue();
    }

    double
    getMolarFactor(void) const
    {
      return getVolume() * fnd::avogadrosNumber;
    }

    // Inserts elements conveying the compartment name, index, and volume.
    void
    insertState(xmlpp::Element* pCompartmentsElt) const
      throw(std::exception);
  };
}

#endif
