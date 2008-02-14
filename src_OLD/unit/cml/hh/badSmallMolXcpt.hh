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

#ifndef CML_BADSMALLMOLXCPT_H
#define CML_BADSMALLMOLXCPT_H

#include "utl/xcpt.hh"
#include "utl/dom.hh"
#include "cml/cptSmallMol.hh"

namespace cml
{
  // Exception thrown when a non-small-mol is called upon to act as a
  // small-mol.
  //
  // This might happen, for instance, if someone tried to do the "nucleotide
  // exchange thing" with a mod-mol.
  class badSmallMolXcpt : 
    public utl::xcpt
  {
    static std::string
    mkMsg(const cptMol* pMol,
	  const xmlpp::Node* pOffendingNode = 0);

  public:
    badSmallMolXcpt(const cptMol* pMol,
		    const xmlpp::Node* pOffendingNode = 0) :
      utl::xcpt(mkMsg(pMol,
		      pOffendingNode))
    {}
  };

  cptSmallMol*
  mustBeSmallMol(cptMol* pMol,
		 const xmlpp::Node* pRequestingNode = 0)
    throw(utl::xcpt);

  const cptSmallMol*
  mustBeSmallMol(const cptMol* pMol,
		 const xmlpp::Node* pRequestingNode = 0)
    throw(utl::xcpt);
}

#endif // CML_BADSMALLMOLXCPT_H
