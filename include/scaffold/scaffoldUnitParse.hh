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

#ifndef SCAFFOLDUNITPARSE_H
#define SCAFFOLDUNITPARSE_H

namespace scaf
{
  // For parsing the reaction generators in the scaffold unit, which
  // about all that the scaffold unit does.

  class addTwentyElevenRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    
  public:
    addTwentyElevenRxnFam(mzr::mzrUnit& refMzrUnit,
			  bnd::molUnit& refMolUnit,
			  plx::plexUnit& refPlexUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}
    
    void
    operator()(xmlpp::Node* pTwentyElevenRxnGenNode) throw(std::exception);
  };

  class addElevenSevenRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    
  public:
    addElevenSevenRxnFam(mzr::mzrUnit& refMzrUnit,
			 bnd::molUnit& refMolUnit,
			 plx::plexUnit& refPlexUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}
    
    void
    operator()(xmlpp::Node* pElevenSevenRxnGenNode) throw(std::exception);
  };

  class addSevenThreeRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    
  public:
    addSevenThreeRxnFam(mzr::mzrUnit& refMzrUnit,
			bnd::molUnit& refMolUnit,
			plx::plexUnit& refPlexUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}
    
    void
    operator()(xmlpp::Node* pSevenThreeRxnGenNode) throw(std::exception);
  };

  class addThreeSevenRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;
    
  public:
    addThreeSevenRxnFam(mzr::mzrUnit& refMzrUnit,
			bnd::molUnit& refMolUnit,
			plx::plexUnit& refPlexUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}
    
    void
    operator()(xmlpp::Node* pThreeSevenRxnGenNode) throw(std::exception);
  };
}

#endif // SCAFFOLDUNITPARSE_H
