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

#ifndef MODKINASEUNITPARSE_H
#define MODKINASEUNITPARSE_H

#include "domUtils/domUtils.hh"
#include "mol/modMol.hh"
#include "stoch/stochUnit.hh"

namespace kinase
{
  // This may require some attention: it's used in various different contexts to
  // generate "activity masks", just setting all the positions in the
  // vector<bool> to true for all the modification sites named in the "name"
  // attribute of the xmlpp::Node*.  Each name must really be the name of
  // a modification site.  No checking on length of mask.
  class setModSiteActive :
    public std::unary_function<xmlpp::Node*, void>
  {
    std::vector<bool>& rMask;
    bnd::modMol* pModMol;
  public:
    setModSiteActive(std::vector<bool>& rActivityMask,
		     bnd::modMol* pModMol) :
      rMask(rActivityMask),
      pModMol(pModMol)
    {}

    void
    operator()(xmlpp::Node* pPhosModSiteRefNode)
      throw(std::exception);
  };

  class addNucBindRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    stoch::stochUnit& rStochUnit;
    
  public:
    addNucBindRxnFam(mzr::mzrUnit& refMzrUnit,
		     bnd::molUnit& refMolUnit,
		     stoch::stochUnit& refStochUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rStochUnit(refStochUnit)
    {}
    
    void
    operator()(xmlpp::Node* pNucBindRxnGenNode) const
      throw(std::exception);
  };

  class addKinaseRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    
  public:
    addKinaseRxnFam(mzr::mzrUnit& refMzrUnit,
		    bnd::molUnit& refMolUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit)
    {}
    
    void
    operator()(xmlpp::Node* pKinaseGenNode) const
      throw(std::exception);
  };

  class addPtaseRxnFam :
    public std::unary_function<xmlpp::Node*, void>
  {
    mzr::mzrUnit& rMzrUnit;
    bnd::molUnit& rMolUnit;
    stoch::stochUnit& rStochUnit;
    
  public:
    addPtaseRxnFam(mzr::mzrUnit& refMzrUnit,
		   bnd::molUnit& refMolUnit,
		   stoch::stochUnit& refStochUnit) :
      rMzrUnit(refMzrUnit),
      rMolUnit(refMolUnit),
      rStochUnit(refStochUnit)
    {}
    
    void
    operator()(xmlpp::Node* pPtaseGenNode) const
      throw(std::exception);
  };
}

#endif // MODKINASEUNITPARSE_H
