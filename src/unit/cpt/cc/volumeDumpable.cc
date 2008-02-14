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

#include "cpt/compartmentGraph.hh"
#include "cpt/volumeDumpable.hh"
#include "cpt/cptEltName.hh"

namespace cpt
{
  volumeDumpable::
  volumeDumpable(const compartmentGraph& rCompartmentGraph) :
    fnd::dumpable<globalDumpArg>(eltName::statStream_volume),
    rGraph(rCompartmentGraph)
  {}

  class addDumpedCompartmentVol :
    public std::unary_function<int, void>
  {
    const compartmentGraph& rGraph;
    double& rTotal;

  public:
    addDumpedCompartmentVol(const compartmentGraph& rCompartmentGraph,
			    double& rTotalVolume) :
      rGraph(rCompartmentGraph),
      rTotal(rTotalVolume)
    {
    }

    void
    operator()(int compartmentIndex) const
    {
      rTotal += rGraph.compartments[compartmentIndex]->getVolume();
    }
  };
    
  void
  volumeDumpable::
  doDump(const globalDumpArg& rDumpArg) const
  {
    std::ostream& rOstream = rDumpArg.getOstream();

    if(rDumpArg.dumpTotal)
      {
	// Emit the total volume of the compartments in the dumparg.
	double totalVol = 0.0;

	std::for_each(rDumpArg.dumpCompartments.end(),
		      rDumpArg.dumpCompartments.end(),
		      addDumpedCompartmentVol(rGraph,
					      totalVol));

	rOstream << totalVol;
      }
    else
      {
	// Emit the volumes of the compartments given in the dumparg.
	std::vector<int>::const_iterator 
	  iCptNdx = rDumpArg.dumpCompartments.begin();

	if(rDumpArg.dumpCompartments.end() != iCptNdx)
	  {
	    double compartmentVol
	      = rGraph.compartments[*iCptNdx]->getVolume();
	      
	    rOstream << compartmentVol;

	    while(rDumpArg.dumpCompartments.end() != ++iCptNdx)
	      {
		double compartmentVol
		  = rGraph.compartments[*iCptNdx]->getVolume();
	      
		rOstream << "\t"
			 <<compartmentVol;
	      }
	  }
      }
  }

  void
  volumeDumpable::
  dumpHeader(const globalDumpArg& rDumpArg) const
  {
    std::ostream& rOstream = rDumpArg.getOstream();

    if(rDumpArg.dumpTotal)
      {
	rOstream << getName();
      }
    else
      {
	std::vector<int>::const_iterator
	  iCptNdx = rDumpArg.dumpCompartments.begin();

	if(rDumpArg.dumpCompartments.end() != iCptNdx)
	  {
	    std::string compartmentName
	      = rGraph.compartments[*iCptNdx]->getName();

	    rOstream << getName()
		     << ":"
		     << compartmentName;

	    while(rDumpArg.dumpCompartments.end() != ++iCptNdx)
	      {
		std::string compartmentName
		  = rGraph.compartments[*iCptNdx]->getName();

		rOstream << "\t"
			 << getName()
			 << ":"
			 << compartmentName;
	      }
	  }
      }
  }
}
