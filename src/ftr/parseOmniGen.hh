/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001, 2008 The Molecular Sciences Institute.
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Original Author:
//   Larry Lok, Research Fellow, Molecular Sciences Institute, 2001

//                     Email: lok@molsci.org
//   
/////////////////////////////////////////////////////////////////////////////

#ifndef FTR_PARSEOMNIGEN_H
#define FTR_PARSEOMNIGEN_H

#include "utl/dom.hh"
#include "mzr/mzrUnit.hh"
#include "mol/molUnit.hh"
#include "plex/plexUnit.hh"

namespace ftr
{
  class parseOmniGen :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    plx::plexUnit& rPlexUnit;

  public:
    parseOmniGen(mzr::mzrUnit& refMzrUnit,
		 bnd::molUnit& refMolUnit,
		 plx::plexUnit& refPlexUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rPlexUnit(refPlexUnit)
    {}

    void
    operator()(xmlpp::Node* pOmniGenNode) const
      throw(utl::xcpt);
  };
}

#endif // FTR_PARSEOMNIGEN_H
