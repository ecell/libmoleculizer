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

#ifndef NUCEXUNITPARSE_H
#define NUCEXUNITPARSE_H

namespace nucEx
{
  // For parsing the reaction generators in the nucleotide exchange unit,
  // which is the main business of this unit.

  class addNucExRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    stoch::stochUnit& rStochUnit;
   
  public:
    addNucExRxnFam(mzr::mzrUnit& refMzrUnit,
		   bnd::molUnit& refMolUnit,
		   plx::plexUnit& refPlexUnit,
		   stoch::stochUnit& refStochUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      rStochUnit(refStochUnit)
    {}
    
    void
    operator()(xmlpp::Node* pNucExRxnGenNode) throw(std::exception);
  };

  class addAutHydRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    stoch::stochUnit& rStochUnit;
    
  public:
    addAutHydRxnFam(mzr::mzrUnit& refMzrUnit,
		    bnd::molUnit& refMolUnit,
		    plx::plexUnit& refPlexUnit,
		    stoch::stochUnit& refStochUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      rStochUnit(refStochUnit)
    {}
    
    void
    operator()(xmlpp::Node* pAutHydGenNode) throw(std::exception);
  };

  class addHetHydRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    stoch::stochUnit& rStochUnit;
    
  public:
    addHetHydRxnFam(mzr::mzrUnit& refMzrUnit,
		    bnd::molUnit& refMolUnit,
		    plx::plexUnit& refPlexUnit,
		    stoch::stochUnit& refStochUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit),
      rStochUnit(refStochUnit)
    {}
    
    void
    operator()(xmlpp::Node* pHetHydGenNode) throw(std::exception);
  };
}

#endif // NUCEXUNITPARSE_H
