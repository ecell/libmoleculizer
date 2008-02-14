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

#ifndef CPT_GLOBALSPECIESDUMPABLEAUX_H
#define CPT_GLOBALSPECIESDUMPABLEAUX_H

namespace cpt
{
  // For adding up the populations of a fixed globalSpecies
  // over several compartments.
  //
  // Used in the templates below and in cptSpeciesDumpable.cc.
  class addDumpedCompartmentPop :
    public std::unary_function<int, void>
  {
    const globalSpecies* pSpecies;
    int& rTotal;

  public:
    addDumpedCompartmentPop(const globalSpecies* pGlobalSpecies,
			    int& rTotalPop) :
      pSpecies(pGlobalSpecies),
      rTotal(rTotalPop)
    {}

    void
    operator()(int compartmentIndex) const
    {
      rTotal += pSpecies->getCompartmentSpecies(compartmentIndex)->getPop();
    }
  };

  // For totalling up the populations of a several globalSpecies
  // over the compartments specified in the dumpArg.
  //
  // Used in the templates below and in cptSpeciesDumpable.cc.
  class addPopsForSpecies :
    public std::unary_function<globalSpecies*, void>
  {
    const globalDumpArg& rDumpArg;
    int& rTotal;

  public:
    addPopsForSpecies(const globalDumpArg& refDumpArg,
		      int& rTotalPop) :
      rDumpArg(refDumpArg),
      rTotal(rTotalPop)
    {}

    void
    operator()(const globalSpecies* pGlobalSpecies) const
    {
      std::for_each(rDumpArg.dumpCompartments.begin(),
		    rDumpArg.dumpCompartments.end(),
		    addDumpedCompartmentPop(pGlobalSpecies,
					    rTotal));
    }
  };

  // For totalling up the populations of several species in a fixed
  // compartment.
  //
  // Used in the templates below and in cptSpeciesDumpable.cc.
  class addSpeciesCompartmentPop :
    public std::unary_function<const globalSpecies*, void>
  {
    int cptNdx;
    int& rTotal;
  public:
    addSpeciesCompartmentPop(int compartmentIndex,
			     int& rTotalPop) :
      cptNdx(compartmentIndex),
      rTotal(rTotalPop)
    {}

    void
    operator()(const globalSpecies* pGlobalSpecies) const
    {
      rTotal += pGlobalSpecies->getCompartmentSpecies(cptNdx)->getPop();
    }
  };
}

#endif // CPT_GLOBALSPECIESDUMPABLEAUX_H
