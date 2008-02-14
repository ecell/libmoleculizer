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

#ifndef CPT_GLOBALDUMPARG_H
#define CPT_GLOBALDUMPARG_H

#include "fnd/basicDumpable.hh"

namespace cpt
{
  class globalDumpArg :
    public fnd::basicDumpable::dumpArg
  {
  public:

    // Compartments whose populations might be used in the dumped figures or
    // headers.
    std::vector<int> dumpCompartments;

    // Whether or not to total the populations in the above compartments,
    // emitting one number (or header) for all of them, instead of one number
    // (or header) for each of them.
    bool dumpTotal;

//     globalDumpArg(std::ostream& rOstream,
// 		  bool doDumpTotal = true) :
//       fnd::basicDumpable::dumpArg(rOstream),
//       dumpTotal(doDumpTotal)
//     {}

    // For conveniently dumping the population of all compartments,
    // either totalled or not.
    globalDumpArg(std::ostream& rOstream,
		  compartmentGraph& rGraph,
		  bool doDumpTotal = true):
      fnd::basicDumpable::dumpArg(rOstream),
      dumpTotal(doDumpTotal)
    {
      // Perhaps I should try to make a habit of using size_t.
      for(int cptNdx = 0;
	  cptNdx < (int) rGraph.compartments.size();
	  ++cptNdx)
	{
	  dumpCompartments.push_back(cptNdx);
	}
    }


    // For dumping the populations of selected compartments,
    // either totalled or not.
    globalDumpArg(std::ostream& rOstream,
		  const std::vector<int>& rDumpCompartments,
		  bool doDumpTotal = true) :
      fnd::basicDumpable::dumpArg(rOstream),
      dumpCompartments(rDumpCompartments),
      dumpTotal(doDumpTotal)
    {}
  };
}

#endif // CPT_GLOBALDUMPARG_H
